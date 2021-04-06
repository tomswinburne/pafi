#include <iostream>
#include <sstream>
#include <iomanip>

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <list>
#include <cmath>

// PAFI files
#include "ConstantsTypes.hpp"
#include "Parser.hpp"
#include "Master.hpp"

#include "StandardSimulator.hpp"
#include "StandardGatherer.hpp"

int main(int narg, char **arg) {
MPI_Init(&narg,&arg);
  MPI_Comm world=MPI_COMM_WORLD;
  run<StandardSimulator,StandardGatherer>(world,"./config.xml");
  MPI_Finalize();
};
