#include "pafi.hpp"

int main(int narg, char **arg)
{
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
  MPI_Group world;
  MPI_Comm_group(MPI_COMM_WORLD,&world);

  if(rank==0) std::cout << "PAFI: MPI_Init done\n";

  // Parse input file
  Parser params("./config.xml");

  if(rank==0) std::cout << "PAFI: XML read\n";

  // Check proposed MPI division is possible
  if(nProcs % params.CoresPerWorker) {
    if(rank == 0) printf("ERROR: nProcs % CoresPerWorker != 0\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int nInstance = nProcs / params.CoresPerWorker;  // workers per iteration

  // Set up MPI communicators for each instance
  int instance = rank / (params.CoresPerWorker);

  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  Simulator sys(instance_comm,params,instance);
  if(rank==0) std::cout<<"SYS LOADED\n";

  sys.make_path(params.KnotList);
  if(rank==0) {
    std::cout<<"PATH MADE\n";
    //std::ofstream out;
    //std::string fn;
    //fn = "PathNormdNorm.dat";
    //out.open(fn.c_str(),std::ofstream::out);
    //out.close();
  }
  std::vector<double> results;
  double T = boost::lexical_cast<double>(params.parameters["LowTemperature"]);
  double position = boost::lexical_cast<double>(params.parameters["position"]);
  //std::vector<std::pair<double,double>> rE;

  // also thermalizes
  std::cout<<std::setprecision(12)<<std::scientific;
  double sT,aT,r=0.,dr = 1./((double)params.nPlanes);

  double E=0.,sE=0.;

  results.clear();
  sys.setup(r,T);
  sT = sys.thermalize();
  aT = sys.sample(results);

  sE = sys.getEnergy();


  while ( r <= 1.+0.01 ) {
    results.clear();
    sys.setup(r,T);
    sT = sys.thermalize();
    aT = sys.sample(results);
    //E += - dr * results[0];
    //if(rank==0) std::cout<<r<<" "<<sys.getEnergy()-sE<<" "<<E<<"\n";
    r += dr;
  }


  //sys.evaluate(position,T,results);


  /*

  //std::vector<double> x(3*sys.natoms,0.);
  //lammps_gather_atoms(sys.lmp,(char *) "norm",1,3,&x[0]);


  //Integrator calculator(sys,params);

  // runs the vector of strings through LAMMPS
  sys.run_vector(params.RunScript);

  double finaltemp = sys.temperature();
  double *temps = new double[nInstance];

  for (int i = 0; i < nInstance; i++) temps[i] = 0.0;
  int instance_rank;

  MPI_Comm_rank(instance_comm,&instance_rank);
  if (instance_rank == 0) temps[instance] = finaltemp;

  double *alltemps = new double[nInstance];
  MPI_Allreduce(temps,alltemps,nInstance,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


  if (rank == 0) for (int i = 0; i < nInstance; i++) \
    printf("Instance %d, final temp = %g\n",i+1,alltemps[i]);

  delete [] temps;
  delete [] alltemps;



  int imax = 0;
  double max = 0.;
  std::vector<double> ox,x,f;
  ox = std::vector<double>(3*natoms,0.);
  x = std::vector<double>(3*natoms,0.);

  lammps_gather_atoms(lmp,(char *) "x",1,3,&ox[0]);
  ox[37+37*instance] += 0.1;
  lammps_scatter_atoms(lmp,(char *) "x",1,3,&ox[0]);
  ox[37+37*instance] -= 0.1;
  */

  // close down LAMMPS instances
  sys.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
