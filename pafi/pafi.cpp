#include "pafi.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Parse input file
  Parser params("./config.xml");

  double T = boost::lexical_cast<double>(params.parameters["LowTemperature"]);


  if(rank==0) std::cout << "PAFI: XML read\n";

  if(nProcs%params.CoresPerWorker!=0) {
    if(rank==0) std::cout<<"CoresPerWorker must factorize nProcs!\n";
    exit(-1);
  }

  const int nWorkers = nProcs / params.CoresPerWorker;
  const int instance = rank / params.CoresPerWorker;
  const int local_rank = rank % params.CoresPerWorker;
  const int nResults = 7;


  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);
  Simulator sim(instance_comm,params,rank);

  sim.make_path(params.KnotList);
  if(rank==0) std::cout<<"\n\nPATH LOADED\n";

  const int vsize = 3 * sim.natoms;


  double *all_results = new double[nResults*nWorkers];
  double *local_results = new double[nResults*nWorkers];
  double *results = new double[nResults];
  std::vector<double> local_deviation(vsize,0.);
  std::vector<double> all_deviation(vsize,0.);
  std::vector<double> all_deviation_sq(vsize,0.);


  for (int i=0;i<nResults;i++) results[i] = 0.;
  for (int i=0;i<nResults*nWorkers;i++) local_results[i] = 0.;
  for (int i=0;i<nResults*nWorkers;i++) all_results[i] = 0.;



  for (double r=0.005; r<=1.; r+=0.1 ) {
    for (int i=0;i<nResults;i++) results[i] = 0.;
    for (int i=0;i<nResults*nWorkers;i++) local_results[i] = 0.;
    for (int i=0;i<nResults*nWorkers;i++) all_results[i] = 0.;

    local_deviation = std::vector<double>(vsize,0.);

    sim.sample(r,T,results, local_deviation);

    if(local_rank == 0) for(int i=0;i<nResults;i++)
      local_results[instance*nResults + i] = results[i];

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(local_results,all_results,nResults*nWorkers,
      MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Gather(&local_deviation[0],vsize,MPI_DOUBLE,&all_deviation[0],vsize,
      MPI_DOUBLE,0,MPI_COMM_WORLD);


    if (rank==0) {
      std::cout<<"SAMPLED: r:"<<r<<" T:"<<T<<"\n";
      for(int i=0;i<nWorkers;i++) {
        for(int j=0;j<nResults;j++) std::cout<<all_results[i*nResults+j]<<" ";
        std::cout<<"\n";
      }
      std::cout<<"\n\n\n";
    }
    //sim.write(r,"PDN_"+boost::lexical_cast<std::string>(r));
  }

  // close down LAMMPS instances
  sim.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
