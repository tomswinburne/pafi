#include "pafi-lammps-preparation.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Load input file
  Parser params("./config.xml",false);
  params.CoresPerWorker = nProcs;

  if(nProcs>params.CoresPerWorker) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"pafi-lammps-path should only be run with a single worker!\n";
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

  int fileindex=1;
  std::string cmd;
  double dr,E_init,nm,E_tot,fs,dF_int=0.0,dF,F_bar=0.0,E_bar=0.0,psi;
  double *f,*lmp_ptr;

  std::vector<double> dFa,ra,dEa;

  params.seed(123);
  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,0,0,&instance_comm);

  Simulator sim(instance_comm,params,instance,nRes);
  if(!sim.has_pafi) {
    if(rank==0)
      std::cout<<"PAFI Error: missing "<<sim.pafi_package<<" package in LAMMPS"<<std::endl;
    exit(-1);
  }
  if(sim.error_count>0 && local_rank==0) std::cout<<sim.last_error()<<std::endl;

  auto cell = sim.getCellData();

  if(rank==0) {
    std::cout<<"Loaded input data of "<<sim.natoms<<" atoms\nSupercell Matrix:\n";
    for(int i=0;i<3;i++) {
      std::cout<<"\t";
      for(int j=0;j<3;j++) std::cout<<sim.pbc.cell[i][j]<<" ";
      std::cout<<"\n";
    }
    std::cout<<"\n\n";
  }

  sim.make_path(params.PathwayConfigurations,params.real_coord);

  std::vector<double> sample_r = params.sample_r(sim.pathway_r);

  f = new double[3*sim.natoms];

  
  if(rank==0) {
    std::cout<<"\n\nPath Loaded\n\n";
    std::cout<<std::setw(17)<<"              r";
    std::cout<<std::setw(17)<<"          index";
    std::cout<<std::setw(19)<<"        E - E_0";
    std::cout<<std::setw(17)<<"         -dF/dr"<<std::endl;
  }

  for (auto r : sample_r) {

    sim.populate(r,nm,0.0);

    sim.run_script("PreRun");  // Stress Fixes

    // pafi fix
    cmd = "run 0\n"; // to ensure the PreRun script is executed
    cmd += "run 0\nfix hp all pafi pafi_path 0.0 ";
    cmd += params.parameters["Friction"]+" ";
    cmd += params.seed_str()+" overdamped 1 com 0\n run 0";
    sim.run_commands(cmd);

    if(params.preMin) {
      #ifdef VERBOSE
      if(rank==0)
        std::cout<<"LAMMPSSimulator.populate(): minimizing"<<std::endl;
      #endif
      cmd = "min_style fire\n minimize 0 0.01 ";
      cmd += params.parameters["MinSteps"]+" "+params.parameters["MinSteps"];
      sim.run_commands(cmd);
    }

    cmd = "run 0";
    sim.run_commands(cmd);

    E_tot = sim.getForceEnergy(f);


    fs=0.0;
    for(int i=0; i<3*sim.natoms; i++) fs += f[i]*f[i];
    dF = 1.0*(sim.get_fix("hp",1,0));
    psi = 1.0*(sim.get_fix("hp",1,2));

    if(r==params.startr) E_init = E_tot;

    ra.push_back(r);
    dFa.push_back(dF*nm);
    dEa.push_back(E_tot-E_init);

    if(rank==0) {
      std::cout<<std::setw(17)<<r<<" ";
      std::cout<<std::setw(17)<<fileindex<<" ";
      std::cout<<std::setw(17)<<std::setprecision(5)<<E_tot-E_init<<" ";
      std::cout<<std::setw(17)<<std::setprecision(5)<<dF*nm<<std::endl;
    }

    sim.lammps_dump_path("dumps/pafipath."+std::to_string(fileindex)+".data",r);
    fileindex++;

    cmd = "unfix hp";
    sim.run_commands(cmd);
    sim.run_script("PostRun");  // Stress Fixes
  }


  // close down LAMMPS instances
  sim.close();


  spline dFspl,dEspl;
  dFspl.set_points(ra,dFa);
  dEspl.set_points(ra,dEa);
  double dF_dense_int = 0.0,F_bar_dense=0.0,E_bar_dense;
  for (double r = params.startr; r <= params.stopr+0.5*dr/30.; r += dr/30. ) {
    dF_dense_int -= dFspl(r)/2.0 * dr/30.;
    F_bar_dense = std::max(F_bar_dense,dF_dense_int);
    dF_dense_int -= dFspl(r)/2.0 * dr/30.;
    E_bar_dense = std::max(E_bar_dense,dEspl(r));
  }


  if(rank==0) {
    std::cout<<"\n\n\tBARRIER FROM ENERGY: "<<std::setprecision(5)<<E_bar_dense;
    std::cout<<"eV, FORCE INTEGRATION:"<<std::setprecision(5)<<F_bar_dense<<"eV";
    std::cout<<"\n\n\tABSOLUTE ERROR: "<<std::setprecision(5)<<std::fabs(F_bar_dense-E_bar_dense)*1000.<<"meV";
    std::cout<<", RELATIVE ERROR: "<<std::setprecision(5)<<(F_bar_dense/E_bar_dense-1.0)*100.<<"%"<<std::endl;
    if(std::fabs(F_bar_dense-E_bar_dense)<=0.005) {
      std::cout<<"\n\tERROR < 5meV - WITHIN LIKELY POTENTIAL ACCURACY, OK FOR SAMPLING!";
      std::cout<<"\n\n\tTO FURTHER REDUCE ERROR, ";
    } else std::cout<<"\n\tERROR > 5meV, PROBABLY TOO HIGH! ";
    std::cout<<"CONSIDER MORE NEB IMAGES AND/OR INCREASING nPlanes IN config.xml\n\n\n"<<std::endl;;
  }

  MPI_Finalize();

}
