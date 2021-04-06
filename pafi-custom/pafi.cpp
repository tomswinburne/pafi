#include "Master.hpp"

#include "CustomSimulator.hpp"
#include "CustomGatherer.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  MPI_Comm world=MPI_COMM_WORLD;
  run<CustomSimulator,CustomGatherer>(world,"./config.xml");
  MPI_Finalize();
};
