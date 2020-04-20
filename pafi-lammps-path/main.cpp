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
      std::cout<<"Could not write to output path / find directory! Exiting!"<<std::endl;
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
    std::cout<<"Loaded input data of "<<sim.natoms<<" atoms\n";
    std::cout<<"Supercell Matrix:\n";
    auto cell = sim.getCellData();
    for(int i=0;i<3;i++) {
      std::cout<<"\t";
      for(int j=0;j<3;j++) std::cout<<sim.pbc.cell[i][j]<<" ";
      std::cout<<"\n";
    }
    std::cout<<"\n\n";
  }

  sim.make_path(params.KnotList);
  MPI_Barrier(MPI_COMM_WORLD);
  int fileindex=1;
  if(instance==0) {
    double *t,*f;
    double dr,E,nm,fE,fs;
    t = new double[3*sim.natoms];
    f = new double[3*sim.natoms];

    if (params.nPlanes>1) dr = (params.stopr-params.startr)/(double)(params.nPlanes-1);
    else dr = 0.1;

    if(rank==0) std::cout<<"\n\nPath Loaded\n\n";

    for (double r = params.startr; r <= params.stopr+0.5*dr; r += dr ) {
      double scale[3] = {1.0,1.0,1.0};
      sim.populate(r,scale,nm);

      std::string cmd = "run 0";
      sim.run_commands(cmd);

      E = sim.getForceEnergy(f);

      fs=0.0;
      for(int i=0; i<3*sim.natoms; i++) {
        t[i] = sim.pathway[i].deriv(1,r);
        fs += f[i]*f[i];
      }

      if(r==params.startr) fE = E;
      std::cout<<std::setprecision(15)<<r<<" "<<fileindex<<" "<<E-fE<<" "<<nm<<" "<<sqrt(fs)<<std::endl;
      //sim.write_dev("path_dpath/path_dpath_f_"+std::to_string(r),r,t,f);

      sim.lammps_path_write("pafipath."+std::to_string(fileindex)+".data",r);
      fileindex++;
    }
  }


  // close down LAMMPS instances
  sim.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);



  MPI_Finalize();

}
