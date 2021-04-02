#ifndef SIM_H
#define SIM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <set>
#include <list>

#include "mpi.h"
#include "ConstantsTypes.hpp"
#include "Parser.hpp"
#include "Boundary.hpp"
#include "Spline.hpp"


class GenericSimulator {

public:

  GenericSimulator(MPI_Comm &instance_comm, Parser &p, Holder &h, int t);

  virtual void load_config(std::string file_string,double *x){};

  virtual void run_script(std::string sn){};

  virtual void run_commands(std::vector<std::string> strv){};

  virtual void setup(double r, double T){};

  virtual double getEnergy(){
    return 0.0;
  };


  virtual double getForceEnergy(double *f){
    return 0.0;
  };

  virtual std::array<double,9> getCellData(){
    std::array<double,9> res;
    res.fill(0.0);
    return res;
  };

  virtual void populate(double r, double T){};

  virtual void rescale_cell(double T){};

  // LAMMPS INDEPENDENT

  void write(std::string fn, double r);

  void write_dev(std::string fn, double r, double *dev);

  void make_path(std::vector<std::string> knot_list);

  double path(int i, double r, int d, double s);

  void evaluate(std::vector<double> &results);

  void expansion(double T,double *newscale);

  virtual void close() {
    delete [] x;
  };

  MPI_Comm *comm;
  Holder results;
  double scale[3],refE,refT,refP, *x;
  int natoms, tag, nknots, out_width, error_count;
  int local_rank, local_size, nlocal, offset;
  MinImage pbc;
  Parser *parser;
  Holder *params;
  std::vector<spline> pathway, splines;
  std::vector<double> pathway_r;
  bool s_flag,has_pafi,spline_path;
  std::string last_error_message;
private:
  /* nothing */
};

#endif
