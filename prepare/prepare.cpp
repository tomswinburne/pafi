#include "prepare.hpp"

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

  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,0,0,&instance_comm);

  Simulator sim(instance_comm,params,rank);

  if(rank==0) std::cout<<"sim LOADED\n";

  sim.make_path(params.KnotList);
  if(rank==0) std::cout<<"PATH MADE\n";
  std::vector<double> results, deviation;
  double T = boost::lexical_cast<double>(params.parameters["LowTemperature"]);
  double r = boost::lexical_cast<double>(params.parameters["position"]);
  sim.sample(r,T,results,deviation);
  if(rank==0) {
    std::cout<<"SAMPLED: r:"<<r<<" T:"<<T<<"\n";
    for(auto rr: results) std::cout<<rr<<" ";
    std::cout<<std::endl;
  }
  r += 0.1;
  sim.sample(r,T,results,deviation);
  if(rank==0) {
    std::cout<<"SAMPLED: r:"<<r<<" T:"<<T<<"\n";
    for(auto rr: results) std::cout<<rr<<" ";
    std::cout<<std::endl;
  }


  // close down LAMMPS instances
  sim.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
