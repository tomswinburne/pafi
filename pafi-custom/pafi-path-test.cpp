#include "Tester.hpp"

#include "CustomSimulator.hpp"
#include "CustomGatherer.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  MPI_Comm world=MPI_COMM_WORLD;
  test<CustomSimulator,CustomGatherer>(world,"./config.xml",false);
  MPI_Finalize();
};
