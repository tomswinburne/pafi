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

  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  if(rank==0) std::cout<<"\n\nInitializing "<<nWorkers<<" workers "
                      "with "<<parser.CoresPerWorker<<" cores per worker\n\n";
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
  double *dev = new double[vsize*(1+(instance==0))];

  // generic - validity
  int *valid = new int[1+nWorkers*(instance==0)];

  // generic - data
  double p_valid,*data=NULL;
  int total_valid, dsize = -1, raw_data_open = 0;
  DataGatherer g;

  // for testing
  std::vector<double> dF,maxjumpr,dE;


  // dump_files
  double T = 0.0;
  dump_suffix = std::to_string(int(T))+"K_"+std::to_string(dump_index);
  dump_file = parser.dump_dir + "/raw_ensemble_output_"+dump_suffix;
  if(rank==0) raw_data_open = g.initialize(parser,dump_file,nWorkers);
  MPI_Bcast(&raw_data_open,1,MPI_INT,0,MPI_COMM_WORLD);
  if(raw_data_open==0) {
    if(rank==0) std::cout<<"Could not open "<<dump_file<<"! EXIT"<<std::endl;
    exit(-1);
  }

  if(rank==0) {
    std::cout<<"\n\n\n*****************************\n\n";
    std::cout<<"T=0K test run, checking pathway for force integration accuracy\n";
    std::cout<<"\n*****************************\n\n";
    sim.screen_output_header(0.0);
  }

  dsize = -1;
  for(auto r: sim.sample_r) {
    total_valid=0;
    // sample
    for(i=0;i<vsize;i++) dev[i] = 0.0;
    sim.sample(r, T, dev);

    // is it valid
    valid[0] = int(sim.results["Valid"]+0.01);
    MPI_Gather(valid,1,MPI_INT,valid+1,1,MPI_INT,0,MPI_COMM_WORLD);

    // reset deviation if invalid
    if(valid[0]==0) for(i=0;i<vsize;i++) dev[i] = 0.0;
    MPI_Reduce(dev,dev+vsize,vsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    // declare data here once, after first simulation for flexibility
    if(dsize<0) {
      dsize = sim.results.size();
      data = new double[dsize];
      if(rank==0) g.prepare(sim.results);
    }

    i=0; for(auto res: sim.results) data[i++] = res.second;
    MPI_Gather(data,dsize,MPI_DOUBLE,g.all_data,dsize,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(rank==0) total_valid += g.ensemble(r,valid+1);
    MPI_Bcast(&total_valid,1,MPI_INT,0,MPI_COMM_WORLD);

    if(rank==0) {
      if(total_valid>0) {
        for(j=0;j<vsize;j++) dev[j+vsize] /= 1.0*total_valid;
        dev_file = "dev_"+std::to_string(r)+"_"+dump_suffix+".dat";
        sim.write_dev(parser.dump_dir+"/"+dev_file,r,dev+vsize);
      }
      sim.fill_results(r,g.ens_data);
      sim.screen_output_line(r);
      g.next(); // wipe ens_data
      dF.push_back(sim.results["aveF"]);
      dE.push_back(sim.results["MinEnergy"]);
      maxjumpr.push_back(sim.results["MaxJump"]);
    }
  }

  // free_energy_profile
  if(rank==0) {
    g.close();
    std::cout<<"\nT=0K test run complete, testing path....\n\n";

    std::cout<<"Absolute Forces and differences between knots: \n";

    dump_file = parser.dump_dir + "/free_energy_profile_"+dump_suffix;
    double F_bar_dense = sim.integrate(dump_file);
    double E_bar_dense=0.;

    double _dF,_ddF,ddF[2] = {0.,0.},_F=0.;
    bool warning=false;
    int i=0;
    for (;i<dF.size()-1;i++) {
      E_bar_dense = std::max(E_bar_dense,dE[i]-dE[0]);
      _ddF = std::fabs(dF[i+1]-dF[i]);
      _dF = -(dF[i]+dF[i+1])/2. * (sim.sample_r[i+1]-sim.sample_r[i]);

      std::cout<<"\n\t Knot "<<i+1<<": r = "<<sim.sample_r[i]<<" |dF| = "<<std::fabs(dF[i])<<" , F = "<<_F<<", E = "<<dE[i]-dE[0]<<"\n";
      std::cout<<"\n\t\t   Max | X_postmin - X | = "<<maxjumpr[i];
      if(maxjumpr[i]>0.02) {
        warning = true;
        std::cout<<" !! This should be zero for a properly discretized minimum energy path!\n"
        " Values greater than ~0.02 should be considered risky;"
        " Perhaps consider more integration points here. WARNING";
      }

      std::cout<<"\n\t\t        | dF_next - dF | = "<<_ddF;
      std::cout<<"\n\t\t Approx |  F_next - F  | = "<<_dF;
      std::cout<<"\n\t\t        |  r_next - r  | = "<<sim.sample_r[i+1]-sim.sample_r[i]<<"\n";

      if(_ddF<0.02 && std::fabs(dF[i])<0.02) {
        warning = true;
        std::cout<<"Very flat segment! "
        "Perhaps remove knot "<<i+1<<" or "<<i+2<<"?\n. FAIL";
      }

      if(_dF>0.1) {
        warning = true;
        std::cout<<"Large free energy change! "
        "Consider more integration points here. FAIL";
      }

      ddF[0] += _ddF/double(dF.size()-1);
      ddF[1] += _ddF*_ddF/double(dF.size()-1);

      _F += -(dF[i]+dF[i+1])/2. * (sim.sample_r[i+1]-sim.sample_r[i]);
      std::cout<<"\n -------- \n";
    }
    ddF[1] -= ddF[0]*ddF[0];
    std::cout<<"\n\t Knot "<<i+1<<": r = "<<sim.sample_r[i]<<" |dF| = "<<std::fabs(dF[i])<<" , F = "<<_F<<"\n";
    std::cout<<"\n\t\t   Max | X_postmin - X | = "<<maxjumpr[i];
    std::cout<<"\n -------- \n";
    std::cout<<"\n\n\tAverage, Std in d|dF|:"<<ddF[0]<<" , "<<sqrt(ddF[1])<<std::endl;



    std::cout<<"\n\n--------------- Integration checks at zero temperature -----\n";
    std::cout<<"\n\tEnergy Barrier: "<<std::setprecision(5)<<E_bar_dense;
    std::cout<<" eV, Force Integration Barrier: "<<std::setprecision(5)<<F_bar_dense<<" eV";
    std::cout<<"\n\n\tAbsolute error of: "<<std::setprecision(5)<<std::fabs(F_bar_dense-E_bar_dense)*1000.<<" meV ";
    std::cout<<" ("<<std::setprecision(5)<<(F_bar_dense/E_bar_dense-1.0)*100.<<"%)"<<std::endl;
    if(std::fabs(F_bar_dense-E_bar_dense)<=0.003) {
      std::cout<<"\n\tError < 3meV - within likely potential accuracy, great for sampling! OK\n";
    } else if(std::fabs(F_bar_dense-E_bar_dense)<0.01) {
      std::cout<<"\n\tAbove target error of ~3meV, but probably OK for sampling CAUTION\n";
    } else {
      warning=true;
      std::cout<<"\n\tError > 10meV  - if using this setup, the "
      "error of "<<F_bar_dense-E_bar_dense<<" eV MUST be propagated to"
      "PAFI results! FAIL\n";
    }
    if(std::fabs(F_bar_dense-E_bar_dense)>=0.001) {
      std::cout<<"\n\tTo reduce this error: ";
      std::cout<<"\n\t More NEB images and/or increasing nPlanes in config.xml";
      std::cout<<"\n\t Try Rediscretize==0/1 to change how images are interpolated\n\n"<<std::endl;
    }

    if(warning) std::cout<<"\nPathway checks failed / warnings generated! See above\n"<<std::endl;
    else std::cout<<"\nPathway checks passed!\n"<<std::endl;

  }

  // close down LAMMPS instances
  sim.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
