#include "Tester.hpp"

#include "LAMMPSSimulator.hpp"
#include "GenericGatherer.hpp"

int main(int narg, char **arg) {
MPI_Init(&narg,&arg);
  MPI_Comm world=MPI_COMM_WORLD;
  test<LAMMPSSimulator,GenericGatherer>(world,"./config.xml",false);
  MPI_Finalize();
};
