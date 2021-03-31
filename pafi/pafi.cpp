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
  const int min_valid = int(params.redo_thresh*nWorkers*params.nRepeats);

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

  // data holders
  const int vsize = 3 * sim.natoms;
  double *dev = new double[vsize*(1+(instance==0))];
  int *valid = new int[1+nWorkers*(instance==0)];
  double *all_data = NULL, *ens_data=NULL, *data=NULL;
  int total_valid,repeats,dsize=-1;

  for(double T = params.lowT; T <= params.highT;) {

    dump_suffix = std::to_string(int(T))+"K_"+std::to_string(dump_index);
    dump_file = params.dump_dir + "/raw_ensemble_output_"+dump_suffix;

    if(rank==0) {
      raw.open(dump_file.c_str(),std::ofstream::out);
      std::cout<<"\nStarting T="<<T<<"K run\n\n";
      std::cout<<"<> == time averages,  av/err over ensemble"<<std::endl;
      sim.screen_output_header();
    }
    for(auto r: sim.sample_r) {
      total_valid=0;
      repeats=0;
      while((repeats<=params.nRepeats) and (total_valid<=min_valid)) {
        repeats++;
        // sample
        for(i=0;i<vsize;i++) dev[i] = 0.0;
        sim.sample(r, T, dev);

        // is it valid
        valid[0] = int(sim.results["Valid"]+0.01);
        MPI_Gather(valid,1,MPI_INT,valid+1,1,MPI_INT,0,MPI_COMM_WORLD);

        // reset deviation if invalid
        if(valid[0]==0) for(i=0;i<vsize;i++) dev[i] = 0.0;
        MPI_Reduce(dev,dev+vsize,vsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

        // declare everything here for flexibility
        if(dsize<0) {
          dsize=sim.results.size();
          data = new double[dsize];
          if(rank==0) {
            i=1;
            raw<<"# 0: r ";
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
          raw<<r<<" ";
          for(i=0;i<dsize*nWorkers;i++) raw<<all_data[i]<<" ";
          // ensemble average
          for(j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
          for(i=0;i<nWorkers;i++) if(valid[1+i]) {
            total_valid++;
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
          sim.fill_results(r,ens_data);
          sim.screen_output_line(r);
        }
      }
    }

    if(rank==0) {
      raw.close();
      dump_file = params.dump_dir + "/free_energy_profile_"+dump_suffix;
      double barrier = sim.integrate(dump_file);
      std::cout<<"T="<<T<<"K run complete; est. barrier: "<<barrier<<"eV"<<std::endl;
    }

    if( params.highT > params.lowT + 1.0 && params.TSteps > 1) \
      T += ( params.highT - params.lowT ) / ( params.TSteps - 1 );
    else break;
  }

  // close down LAMMPS instances
  sim.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
