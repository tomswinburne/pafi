#include "pafi.hpp"

int main(int narg, char **arg) {

  MPI_Init(&narg,&arg);
  int rank, nProcs, i, j;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // ************************ READ CONFIG FILE **********************************
  Parser params("./config.xml",false);

  if(!params.xml_success) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Configuration file could not be read!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }

  if(nProcs%params.CoresPerWorker!=0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"CoresPerWorker must factorize nProcs!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }
  // ************************ READ CONFIG FILE ***********************************


  // ************************ DUMP FOLDER *********************************
  std::ofstream raw;
  int dump_index=-1;
  std::string dump_suffix, dump_file, dev_file;
  if(rank==0) params.find_dump_file(raw,dump_index);
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
  // ************************ DUMP FOLDER *********************************



  // ******************* SET UP WORKERS ***************************************
  MPI_Barrier(MPI_COMM_WORLD);
  const int nWorkers = nProcs / params.CoresPerWorker;
  const int instance = rank / params.CoresPerWorker;
  const int local_rank = rank % params.CoresPerWorker;

  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  if(rank==0) std::cout<<"\n\nInitializing "<<nWorkers<<" workers "
                      "with "<<params.CoresPerWorker<<" cores per worker\n\n";
  // see GlobalSeed
  params.seed(instance);

  Simulator sim(instance_comm,params,instance);
  if(!sim.has_pafi) exit(-1);
  if(rank==0)  std::cout<<params.welcome_message();

  sim.make_path(params.PathwayConfigurations);
  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";
  // ******************* SET UP WORKERS ****************************************


  // ********************** SAMPLING *******************************************

  const int vsize = 3 * sim.natoms;
  double *dev = new double[vsize*(1+(instance==0))];
  int *valid = new int[1+nWorkers*(instance==0)];
  double *all_data = NULL, *ens_data=NULL, *data=NULL;
  int dsize=-1;
  std::vector<double> dF,maxjumpr,dE;

  double T = 0.0;
  dump_suffix = std::to_string(int(T))+"K_"+std::to_string(dump_index);
  dump_file = params.dump_dir + "/raw_ensemble_output_"+dump_suffix;

  if(rank==0) {
    raw.open(dump_file.c_str(),std::ofstream::out);
    std::cout<<"\n\n\n*****************************\n\n";
    std::cout<<"T=0K test run, checking pathway for force integration accuracy\n";
    std::cout<<"\n*****************************\n\n";
    sim.screen_output_header();
    dF.clear();
    maxjumpr.clear();
    dE.clear();
  }
  for(auto r: sim.sample_r) {
    // sample
    for(i=0;i<vsize;i++) dev[i] = 0.0;
    sim.sample(r, T, dev);

    // is it valid
    valid[0] = int(sim.results["Valid"]+0.01);
    MPI_Gather(valid,1,MPI_INT,valid+1,1,MPI_INT,0,MPI_COMM_WORLD);

    // reset deviation if invalid
    if(valid[instance]==0) for(i=0;i<vsize;i++) dev[i] = 0.0;
    MPI_Reduce(dev,dev+vsize,vsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    // declare everything here for flexibility
    if(dsize<0) {
      dsize=sim.results.size();
      data = new double[dsize];
      if(rank==0) {
        i=0;
        raw<<"# ";
        for(auto res: sim.results) raw<<i++<<": "<<res.first<<"  ";
        raw<<std::endl;
        all_data = new double[dsize*nWorkers];
        ens_data = new double[dsize*2+1];
      }
    }
    i=0; for(auto res: sim.results) data[i++] = res.second;
    MPI_Gather(data,dsize,MPI_DOUBLE,all_data,dsize,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(rank==0) {
      // raw output
      for(i=0;i<dsize*nWorkers;i++) raw<<all_data[i]<<" ";
      raw<<std::endl;
      // ensemble average
      for(j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
      for(i=0;i<nWorkers;i++) if(valid[1+i]) {
        ens_data[2*dsize] += 1.0;
        for(j=0;j<dsize;j++) {
          ens_data[j]+=all_data[i*dsize+j];
          ens_data[j+dsize] += all_data[i*dsize+j] * all_data[i*dsize+j];
        }
      }
      if(ens_data[2*dsize]>0.5) {
        for(j=0;j<2*dsize;j++) ens_data[j] /= ens_data[2*dsize];
        for(j=0;j<dsize;j++) ens_data[j+dsize] -= ens_data[j] * ens_data[j];
        // N^2 for average-of-averages
        for(j=0;j<vsize;j++) dev[j+vsize] /= ens_data[2*dsize] * ens_data[2*dsize];
        dev_file = "dev_"+std::to_string(r)+"_"+dump_suffix+".dat";
        sim.write_dev(params.dump_dir+"/"+dev_file,r,dev+vsize);
      }
      sim.fill_results(ens_data);
      sim.screen_output_line(r);
      dF.push_back(sim.results["aveF"]);
      dE.push_back(sim.results["MinEnergy"]);
      maxjumpr.push_back(sim.results["MaxJump"]);
    }
  }


  // free_energy_profile
  if(rank==0) {
    raw.close();
    std::cout<<"\nT=0K test run complete, testing path....\n\n";

    std::cout<<"Absolute Forces and differences between knots: \n";

    dump_file = params.dump_dir + "/free_energy_profile_"+dump_suffix;
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
