#include "pafi.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Load input file
  Parser params("./config.xml");


  if(nProcs%params.CoresPerWorker!=0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"CoresPerWorker must factorize nProcs!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }



  // Find fresh dump folder name - no nice solution here as
  // directory creation requires platform dependent features
  // which we omit for portability
  int *int_dump_suffix = new int[1];
  std::ofstream raw;
  std::string params_file = params.dump_dir+"/params_"+std::to_string(int_dump_suffix[0]);

  // try to write to a file with a unique suffix
  if(rank==0) {
    for (int_dump_suffix[0]=0; int_dump_suffix[0] < 100; int_dump_suffix[0]++) {
      params_file = params.dump_dir+"/params_"+std::to_string(int_dump_suffix[0]);
      if(!file_exists(params_file)) {
        raw.open(params_file.c_str(),std::ofstream::out);
        if(raw.is_open()) {
          raw<<params.welcome_message();
          raw.close();
          break;
        }
      }
    }
  }
  MPI_Bcast(int_dump_suffix,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if(int_dump_suffix[0]==100) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Could not find/write to output path \""<<params.dump_dir<<"\" Exiting!"<<std::endl;
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }

  if(rank==0)  std::cout<<params.welcome_message();

  MPI_Barrier(MPI_COMM_WORLD);

  const int nWorkers = nProcs / params.CoresPerWorker;
  const int instance = rank / params.CoresPerWorker;
  const int local_rank = rank % params.CoresPerWorker;
  const int nRes = 7; // TODO remove this
	const int nRepeats = params.nRepeats;



  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  int *seed = new int[nWorkers];
  if(rank==0) for(int i=0;i<nWorkers;i++) seed[i]= 137*i + static_cast<int>(std::time(0))/1000;
  MPI_Bcast(seed,nWorkers,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  params.seed(seed[instance]);

  if(rank==0) std::cout<<"\n\nSet up "<<nWorkers<<" workers with   "<<params.CoresPerWorker<<" cores per worker\n\n";

  Simulator sim(instance_comm,params,rank);

  if (rank == 0) {
    std::cout<<"Loaded input data of "<<sim.getNatoms()<<" atoms\n";
    std::cout<<"Supercell Matrix:\n";
    auto cell = sim.getCellData();
    std::cout<<"\t"<<cell[0]<<" "<<cell[3]<<" "<<cell[4]<<"\n\t0 "<<cell[1]<<" "<<cell[5]<<"\n\t0 0 "<<cell[2]<<"\n\n\n";
  }

  sim.make_path(params.KnotList);
  MPI_Barrier(MPI_COMM_WORLD);

  if(instance==0) {
    double *t,*f;
    double dr,E,nm,fE,fs;
    double scale[3] = {1.0,1.0,1.0};

    if (params.nPlanes>1) dr = (params.stopr-params.startr)/(double)(params.nPlanes-1);
    else dr = 0.1;

    if(rank==0) std::cout<<"\n\nPath Loaded\n\n";

    for (double r = params.startr; r < params.stopr+0.5*dr; r += dr ) {

      sim.populate(r,scale,nm);

      E = sim.getEnergy();

      if(r==params.startr) fE = E;
      if(rank==0) std::cout<<std::setprecision(15)<<r<<" "<<E-fE<<" "<<nm<<std::endl;
    }
  }

  // close down LAMMPS instances
  sim.close();

  delete [] seed;
  delete [] int_dump_suffix;

  // close down MPI
  MPI_Comm_free(&instance_comm);



  MPI_Finalize();

}
