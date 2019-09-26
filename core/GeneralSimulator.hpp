#ifndef SIM_H
#define SIM_H

#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <map>

#include "Constants.hpp"
#include "Parser.hpp"
#include "Boundary.hpp"
#include "Spline.hpp"

class GeneralSimulator {

public:

  virtual void load_config(std::string file_string,std::vector<double> &x){};

  virtual void run_script(std::string sn){};

  virtual void run_commands(std::vector<std::string> strv){};

  virtual void setup(double r, double T){};

  virtual void sample(double r, double T, double *results, double *dev){};

  virtual double getEnergy(){};

  virtual std::array<double,9> getCellData(){};

  virtual void close(){};

  virtual void populate(double r){};

  virtual void rescale_cell(double scale){};

  virtual void sample(double r, double T, std::vector<double> &results){};

  // LAMMPS INDEPENDENT

  void write(std::string fn, double r);

  void write_dev(std::string fn, double r, double *dev, double *dev_sq);

  double expansion(double T);

  void make_path(std::vector<std::string> knot_list);

  void evaluate(std::vector<double> &results);

  double refE;
  int natoms, tag, nknots;
  MinImage pbc;
  Parser *params;
  std::vector<spline> pathway;
  bool s_flag;
private:
  /* nothing */
};


#endif
