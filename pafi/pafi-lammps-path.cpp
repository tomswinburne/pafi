#include "pafi.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Load input file
  Parser params("./config.xml");
  params.CoresPerWorker = nProcs;

  if(nProcs>1) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"pafi-lammps-path should only be run in serial!\n";
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

  //if(rank==0)  std::cout<<params.welcome_message();


  const int nWorkers = nProcs / params.CoresPerWorker;
  const int instance = rank / params.CoresPerWorker;
  const int local_rank = rank % params.CoresPerWorker;
  const int nRes = 7; // TODO remove this
	const int nRepeats = params.nRepeats;



  params.seed(123);
  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,0,0,&instance_comm);

  Simulator sim(instance_comm,params,rank);


  std::cout<<"Loaded input data of "<<sim.natoms<<" atoms\n";
  std::cout<<"Supercell Matrix:\n";
  auto cell = sim.getCellData();
  for(int i=0;i<3;i++) {
    std::cout<<"\t";
    for(int j=0;j<3;j++) std::cout<<sim.pbc.cell[i][j]<<" ";
    std::cout<<"\n";
  }
  std::cout<<"\n\n";

  sim.make_path(params.KnotList);
  int fileindex=1;

  double *t,*f;
  double dr,E,nm,fE,fs;
  t = new double[3*sim.natoms];
  f = new double[3*sim.natoms];

  if (params.nPlanes>1) dr = (params.stopr-params.startr)/(double)(params.nPlanes-1);
  else dr = 0.1;

  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";

  for (double r = params.startr; r <= params.stopr+0.5*dr; r += dr ) {

    sim.populate(r,nm);

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

    sim.lammps_dump_path("dumps/pafipath."+std::to_string(fileindex)+".data",r);
    fileindex++;
  }


  // close down LAMMPS instances
  sim.close();

  MPI_Finalize();

}
