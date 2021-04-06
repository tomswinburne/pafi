#include "Master.hpp"
#include "LAMMPSSimulator.hpp"
#include "GenericGatherer.hpp"


int main(int narg, char **arg) {
MPI_Init(&narg,&arg);
  MPI_Comm world=MPI_COMM_WORLD;
  run<LAMMPSSimulator,GenericGatherer>(world,"./config.xml");
  MPI_Finalize();
};
