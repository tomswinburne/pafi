#include "pafi.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Load input file
  Parser params("./config.xml",false);
  params.CoresPerWorker = nProcs;


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

  Simulator sim(instance_comm,params,instance);
  if(!sim.has_pafi) {
    if(rank==0)
      std::cout<<"PAFI Error: missing USER-MISC package in LAMMPS"<<std::endl;
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

  sim.make_path(params.PathwayConfigurations);

  f = new double[3*sim.natoms];

  if (params.nPlanes>1) dr = (params.stopr-params.startr)/(double)(params.nPlanes-1);
  else dr = 0.1;
  std::vector<double> sample_r;
  if(params.spline_path and not params.match_planes) {
    for (double r = params.startr; r <= params.stopr+0.5*dr; r += dr )
      sample_r.push_back(r);
  } else for(auto r: sim.pathway_r) if(r>=0.0 && r<=1.0) sample_r.push_back(r);



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
    cmd += "run 0\nfix hp all pafi __pafipath 0.0 ";
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

    std::cout<<"\n\n******************************************************************\n\n"
    "\tONLY SIMPLE TESTS PERFORMED HERE!\n"
    "\tPLEASE USE pafi-path-test TO TEST PATHWAY FOR FORCE INTEGRATION\n"
    "\n******************************************************************\n\n"<<std::endl;

    std::cout<<"\n\n\tBARRIER FROM ENERGY: "<<std::setprecision(5)<<E_bar_dense;
    std::cout<<"eV, FORCE INTEGRATION:"<<std::setprecision(5)<<F_bar_dense<<"eV";
    std::cout<<"\n\n\tABSOLUTE ERROR: "<<std::setprecision(5)<<std::fabs(F_bar_dense-E_bar_dense)*1000.<<"meV";
    std::cout<<", RELATIVE ERROR: "<<std::setprecision(5)<<(F_bar_dense/E_bar_dense-1.0)*100.<<"%"<<std::endl;


    if(std::fabs(F_bar_dense-E_bar_dense)<=0.003) {
      std::cout<<"\n\tError < 3meV - verify with pafi-path-test\n";
    } else if(std::fabs(F_bar_dense-E_bar_dense)<0.01) {
      std::cout<<"\n\tAbove target error of ~3meV, verify with pafi-path-test\n";
    } else {
      std::cout<<"\n\tError > 10meV  - if using this setup, the "
      "error of "<<F_bar_dense-E_bar_dense<<" eV MUST be propagated to"
      "PAFI results! verify with pafi-path-test\n";
    }
    if(std::fabs(F_bar_dense-E_bar_dense)>=0.001) {
      std::cout<<"\n\tTo reduce this error: ";
      std::cout<<"\n\t More NEB images and/or increasing nPlanes in config.xml";
      std::cout<<"\n\t Try Rediscretize==0/1 to change how images are interpolated\n\n"<<std::endl;
    }

    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;


  }

  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
