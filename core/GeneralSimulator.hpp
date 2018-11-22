#ifndef SIM_H
#define SIM_H

#include <mpi.h>
#include <fstream>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>

#include "Parser.hpp"
#include "Boundary.hpp"

#include "Spline.hpp"

class GeneralSimulator {

public:

  virtual void load_config(std::string file_string,std::vector<double> &x){};

  virtual void run_script(std::string sn){};

  virtual void run_commands(std::vector<std::string> strv){};

  virtual void setup(double r, double T){};

  virtual double thermalize(){};

  virtual double sample(std::vector<double> &results){};

  virtual double getEnergy(){};

  virtual void close(){};

  virtual void populate(double r){};

  virtual void rescale_cell(double scale){};

  virtual void sample(double r, double T, std::vector<double> &results){};

  // LAMMPS INDEPENDENT

  void write(double r,std::string fn);

  double expansion(double T);

  void make_path(std::vector<std::string> knot_list);

  void evaluate(std::vector<double> &results);

  double refE;
  int natoms, tag, nknots;
  MinImage pbc;
  Parser *params;
  std::vector<spline> pathway;
private:
  /* nothing */
};


#endif
