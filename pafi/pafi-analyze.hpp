#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <cmath>
#include <ctime>



// PAFI files
#include "Constants.hpp"
#include "Parser.hpp"
#include "Spline.hpp"

bool file_exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}



#include "LAMMPSSimulator.hpp"
typedef LAMMPSSimulator Simulator;
