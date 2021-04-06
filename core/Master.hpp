
template <class SimulatorTemplate,class GathererTemplate>
void run(MPI_Comm &world,std::string parser_file) {

  int rank, nProcs, i, j;

  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&nProcs);

  // ************************ READ CONFIG FILE **********************************
  Parser parser(parser_file,false);

  if(!parser.xml_success) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Configuration file could not be read!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }

  if(nProcs%parser.CoresPerWorker!=0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"CoresPerWorker must factorize nProcs!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }
  // ************************ READ CONFIG FILE ***********************************


  // ************************ DUMP FOLDER *********************************
  int dump_index=-1;
  if(rank==0) parser.find_dump_file(dump_index);
  MPI_Bcast(&dump_index,1,MPI_INT,0,world);
  MPI_Barrier(world);
  if(dump_index<0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"PAFI could not write to output path : does the directory ";
      std::cout<<parser.dump_dir<<" exist? Exiting!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }
  // ************************ DUMP FOLDER *********************************



  // ******************* SET UP WORKERS ***************************************
  const int nWorkers = nProcs / parser.CoresPerWorker;
  const int instance = rank / parser.CoresPerWorker;
  const int local_rank = rank % parser.CoresPerWorker;
  const int min_valid = int(parser.redo_thresh*nWorkers*parser.nRepeats);

  // LAMMPS communicators
  MPI_Comm instance_comm;
  MPI_Group world_group, ensemble_group;
  MPI_Comm ensemble_comm;

  MPI_Comm_split(world,instance,0,&instance_comm);

  // local_rank==0 communicators - so verbose!
  int *masters = new int[nWorkers];
  for(i=0;i<nWorkers;i++) masters[i] = parser.CoresPerWorker*i;
  const int *cmasters = masters;
  MPI_Comm_group(world, &world_group);
  MPI_Group_incl(world_group, nWorkers, cmasters, &ensemble_group);
  MPI_Comm_create_group(world, ensemble_group, 0, &ensemble_comm);



  // see GlobalSeed
  parser.seed(instance);

  // set up data gatherer
  GathererTemplate g(parser,nWorkers,dump_index,rank);
  MPI_Bcast(&(g.initialized),1,MPI_INT,0,world);
  if(g.initialized==0) exit(-1);



  SimulatorTemplate sim(instance_comm,parser,g.params,instance);
  if(!sim.has_pafi) exit(-1);
  if(rank==0)  std::cout<<parser.welcome_message();

  sim.make_path(parser.PathwayConfigurations);
  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";
  g.special_r_overwrite(sim.pathway_r);

  // ******************* SET UP WORKERS ****************************************


  // ********************** SAMPLING *******************************************

  // generic - deviation
  const int vsize = 3 * sim.natoms;
  double *dev = new double[vsize*(1+(rank==0))];

  // generic - validity
  int *valid = new int[1+nWorkers*(rank==0)];


  int total_valid;

  if(rank==0) {
    std::cout<<"\n\nInitialized "<<nWorkers<<" workers "
                      "with "<<parser.CoresPerWorker<<" cores per worker\n\n";
    std::cout<<"<> == time averages,  av/err over ensemble"<<std::endl;
  }


  g.screen_output_header();

  while(true) {
    total_valid=0;
    for(int repeat=1;repeat<=parser.nRepeats+parser.maxExtraRepeats;repeat++) {
      // sample
      for(i=0;i<vsize;i++) dev[i] = 0.0;

      // need to distribute parameters....?
      sim.sample(g.params, dev); // sim*(parser, dev)

      // is it valid

      valid[0] = int(sim.results["Valid"]+0.01);
      if(MPI_COMM_NULL != ensemble_comm)
        MPI_Gather(valid,1,MPI_INT,valid+1,1,MPI_INT,0,ensemble_comm);

      // reset deviation if invalid
      if(valid[0]==0) for(i=0;i<vsize;i++) dev[i] = 0.0;
      if(MPI_COMM_NULL != ensemble_comm)
        MPI_Reduce(dev,dev+vsize,vsize,MPI_DOUBLE,MPI_SUM,0,ensemble_comm);

      g.prepare(sim.results); // allocate memory if not already done
      if(MPI_COMM_NULL != ensemble_comm)
        MPI_Gather(g.data,g.dsize,MPI_DOUBLE,g.all_data,g.dsize,MPI_DOUBLE,0,ensemble_comm);
      total_valid += g.collate(valid+1);
      MPI_Bcast(&total_valid,1,MPI_INT,0,world);

      if(rank==0 and parser.nRepeats>1) {
          std::cout<<"Repeat "<<repeat<<"/"<<parser.nRepeats<<", Validity ";
          std::cout<<std::setprecision(4);
          std::cout<<100.0/double(nWorkers*repeat) * total_valid<<"% ";
          if((repeat>=parser.nRepeats) and (total_valid<min_valid))
            std::cout<<"Low validity! Running additional repeat!";
          std::cout<<std::endl;
      }
      if((repeat>=parser.nRepeats) and (total_valid>=min_valid)) break;
    }
    if(rank==0 and total_valid>0 and parser.write_dev) {
      for(j=0;j<vsize;j++) dev[j+vsize] /= 1.0*total_valid;
      sim.write_dev(g.dev_file,g.params["ReactionCoordinate"],dev+vsize);
    }
    g.next(); // wipe ens_data
    MPI_Barrier(world);

    if(g.finished()) break;
  }
  if(rank==0) g.close();

  // close down LAMMPS instances
  sim.close();

  MPI_Barrier(world);

  // close down MPI
  MPI_Comm_free(&instance_comm);
  if(MPI_COMM_NULL != ensemble_comm) MPI_Comm_free(&ensemble_comm);
  MPI_Group_free(&world_group);
  MPI_Group_free(&ensemble_group);
};
