#include <iostream>
#include <fstream>

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <regex>
#include <cmath>

// these are LAMMPS include files
#include "lammps/lammps.h"
#include "lammps/input.h"
#include "lammps/atom.h"
#include "lammps/library.h"

// boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random.hpp>
#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>


/*
char *lammps_run_vector(void *ptr, std::vector<std::string> strv){
  LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *) ptr;
  for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
}
*/

// PAFI files
#include "pbc.hpp"
#include "parser.hpp"
#include "spline.hpp"
#include "Simulator.hpp"
#include "integrator.hpp"

typedef LAMMPSSimulator Simulator;
