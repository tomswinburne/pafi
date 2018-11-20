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
#include <cmath>

// these are LAMMPS include files
#include "lammps/lammps.h"
//#include "lammps/input.h"
//#include "lammps/atom.h"
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
//#include <boost/timer/timer.hpp>
#include <boost/foreach.hpp>

#include <Eigen/Dense> // For Supercell and Hessian


// PAFI files
#include "Spline.hpp"

#include "Boundary.hpp"
#include "Parser.hpp"

#include "LAMMPSSimulator.hpp"
typedef LAMMPSSimulator Simulator;
