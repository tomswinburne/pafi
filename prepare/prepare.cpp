#include "prepare.hpp"

int main(int narg, char **arg)
{
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);


  if(rank==0) std::cout << "PAFI: MPI_Init done\n";

  // Parse input file
  Parser params("./config.xml");

  if(rank==0) std::cout << "PAFI: XML read\n";


  if(nProcs%params.CoresPerWorker!=0) {
    if(rank==0) std::cout<<"CoresPerWorker must factorize nProcs!\n";
    exit(-1);
  }

  //MPI_Group world;
  //MPI_Comm_group(MPI_COMM_WORLD,&world);

  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,me/params.CoresPerWorker,0,&instance_comm);

  Simulator sim(instance_comm,params,rank);

  if(rank==0) std::cout<<"sim LOADED\n";

  sim.make_path(params.KnotList);

  std::vector<double> results, deviation;
  double T = boost::lexical_cast<double>(params.parameters["LowTemperature"]);
  //double r = boost::lexical_cast<double>(params.parameters["position"]);

  for (double r=0.0; r<=1.0; r+=0.1 ) {
    sim.sample(r,T,results,deviation);
    if (rank==0) {
      sim.write(r,"PDN_"+boost::lexical_cast<std::string>(r));
      std::cout<<"SAMPLED: r:"<<r<<" T:"<<T<<"\n";
      for(auto rr: results) std::cout<<rr<<" ";
      std::cout<<std::endl;
    }
  }

  // r += 0.1;
  // sim.sample(r,T,results,deviation);
  // if(rank==0) {
  //   std::cout<<"SAMPLED: r:"<<r<<" T:"<<T<<"\n";
  //   for(auto rr: results) std::cout<<rr<<" ";
  //   std::cout<<std::endl;
  // }
  // close down LAMMPS instances
  sim.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
