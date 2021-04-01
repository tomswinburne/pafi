#include "pafi.hpp"

int main(int narg, char **arg) {

  MPI_Init(&narg,&arg);
  int rank, nProcs, i, j;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // ************************ READ CONFIG FILE **********************************
  Parser parser("./config.xml",false);

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
  MPI_Bcast(&dump_index,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(dump_index<0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Could not write to output path / find directory! Exiting!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }
  std::string dump_suffix, dump_file, dev_file;
  // ************************ DUMP FOLDER *********************************



  // ******************* SET UP WORKERS ***************************************
  MPI_Barrier(MPI_COMM_WORLD);
  const int nWorkers = nProcs / parser.CoresPerWorker;
  const int instance = rank / parser.CoresPerWorker;
  const int local_rank = rank % parser.CoresPerWorker;
  const int min_valid = int(parser.redo_thresh*nWorkers*parser.nRepeats);

  // LAMMPS communicators
  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  // local_rank==0 communicators - so verbose!
  int *masters = new int[nWorkers];
  for(i=0;i<nWorkers;i++) masters[i] = parser.CoresPerWorker*i;
  const int *cmasters = masters;
  MPI_Group world_group, ensemble_group;
  MPI_Comm ensemble_comm;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_incl(world_group, nWorkers, cmasters, &ensemble_group);
  MPI_Comm_create_group(MPI_COMM_WORLD, ensemble_group, 0, &ensemble_comm);

  // see GlobalSeed
  parser.seed(instance);

  Simulator sim(instance_comm,parser,instance);
  if(!sim.has_pafi) exit(-1);
  if(rank==0)  std::cout<<parser.welcome_message();

  sim.make_path(parser.PathwayConfigurations);
  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";
  // ******************* SET UP WORKERS ****************************************


  // ********************** SAMPLING *******************************************

  // generic - deviation
  const int vsize = 3 * sim.natoms;
  double *dev = new double[vsize*(1+(rank==0))];

  // generic - validity
  int *valid = new int[1+nWorkers*(rank==0)];

  // generic - data
  double p_valid,*data=NULL,*all_data=NULL;
  int total_valid, dsize = -1, raw_data_open = 0;
  DataGatherer g;

  if(rank==0) std::cout<<"\n\nInitialized "<<nWorkers<<" workers "
                      "with "<<parser.CoresPerWorker<<" cores per worker\n\n";

  // we want generic structure to go through paramater space, with form
  /*
  Standard PAFI has two parameters: r, T
  could extend to stress for example

  Sweep [ T, r, ]....

  generates dump_suffix, raw_file, dump_file
  */

  for(double T = parser.lowT; T <= parser.highT;) {

    // DUMP
    dump_suffix = std::to_string(int(T))+"K_"+std::to_string(dump_index);
    dump_file = parser.dump_dir + "/raw_data_output_"+dump_suffix;

    if(rank==0) raw_data_open = g.initialize(parser,dump_file,nWorkers);
    MPI_Bcast(&raw_data_open,1,MPI_INT,0,MPI_COMM_WORLD);
    if(raw_data_open==0) {
      if(rank==0) std::cout<<"Could not open "<<dump_file<<"! EXIT"<<std::endl;
      exit(-1);
    }
    // DUMP

    if(rank==0) sim.screen_output_header(T);

    dsize = -1;
    for(auto r: sim.sample_r) {
      total_valid=0;
      for(int repeat=1;repeat<=parser.nRepeats+parser.maxExtraRepeats;repeat++){
        // sample
        for(i=0;i<vsize;i++) dev[i] = 0.0;
        DataVec params;
        params["ReactionCoordinate"] = r;
        params["Temperature"] = T;
        sim.sample(params, dev); // sim*(parser, dev)

        // is it valid
        valid[0] = int(sim.results["Valid"]+0.01);
        if(MPI_COMM_NULL != ensemble_comm)
          MPI_Gather(valid,1,MPI_INT,valid+1,1,MPI_INT,0,ensemble_comm);

        // reset deviation if invalid
        if(valid[0]==0) for(i=0;i<vsize;i++) dev[i] = 0.0;
        if(MPI_COMM_NULL != ensemble_comm)
          MPI_Reduce(dev,dev+vsize,vsize,MPI_DOUBLE,MPI_SUM,0,ensemble_comm);

        // declare data here once, after first simulation for flexibility
        if(dsize<0) {
          dsize = sim.results.size();
          data = new double[dsize];
          all_data = new double[dsize*nWorkers];
          if(rank==0) g.prepare(sim.results);
        }

        i=0; for(auto res: sim.results) data[i++] = res.second;
        if(MPI_COMM_NULL != ensemble_comm)
          MPI_Gather(data,dsize,MPI_DOUBLE,all_data,dsize,MPI_DOUBLE,0,ensemble_comm);

        if(rank==0) total_valid += g.ensemble(r,valid+1,all_data);

        MPI_Bcast(&total_valid,1,MPI_INT,0,MPI_COMM_WORLD);
        p_valid = 100.0 / double(nWorkers*repeat) * total_valid;
        if(rank==0 and parser.nRepeats>1) {
            std::cout<<"Repeat "<<repeat<<"/"<<parser.nRepeats<<", Validity ";
            std::cout<<std::setprecision(4)<<p_valid<<"% ";
            if((repeat>=parser.nRepeats) and (total_valid<min_valid))
              std::cout<<"Low validity! Running additional repeat!";
            std::cout<<std::endl;
        }
        if((repeat>=parser.nRepeats) and (total_valid>=min_valid)) break;

      }
      if(rank==0) {
        if(total_valid>0) for(j=0;j<vsize;j++) dev[j+vsize] /= 1.0*total_valid;

        dev_file = "dev_"+std::to_string(r)+"_"+dump_suffix+".dat";
        sim.write_dev(parser.dump_dir+"/"+dev_file,r,dev+vsize);

        sim.fill_results(r,g.ens_data);
        sim.screen_output_line(r);
        g.next(); // wipe ens_data
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    if(rank==0) {
      g.close();
      // to be replaced...
      dump_file = parser.dump_dir + "/free_energy_profile_"+dump_suffix;
      sim.end_of_cycle(dump_file);
    }

    if( parser.highT > parser.lowT + 1.0 && parser.TSteps > 1) \
      T += ( parser.highT - parser.lowT ) / ( parser.TSteps - 1 );
    else break;
  }

  // close down LAMMPS instances
  sim.close();

  MPI_Barrier(MPI_COMM_WORLD);

  // close down MPI
  MPI_Comm_free(&instance_comm);
  if(MPI_COMM_NULL != ensemble_comm) MPI_Comm_free(&ensemble_comm);
  MPI_Group_free(&world_group);
  MPI_Group_free(&ensemble_group);


  MPI_Finalize();

}
