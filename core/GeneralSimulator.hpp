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

  virtual void load_config(std::string file_string,double *x){};

  virtual void run_script(std::string sn){};

  virtual void run_commands(std::vector<std::string> strv){};

  virtual void setup(double r, double T){};

  virtual void sample(double r, double T, double *results, double *dev){};

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

  virtual void close(){};

  virtual void populate(double r, double T){};

  virtual void sample(double r, double T, std::map<std::string,double> &results){};

  virtual void rescale_cell(double T){};


  // LAMMPS INDEPENDENT

  void write(std::string fn, double r);

  void write_dev(std::string fn, double r, double *dev, double *dev_sq);

  void make_path(std::vector<std::string> knot_list, int real_coord);

  double path(int i, double r, int d, double s);

  void evaluate(std::vector<double> &results);

  void expansion(double T,double *newscale);

  double scale[3];
  int natoms, tag, nknots, nres;
  double MinEnergy;
  MinImage pbc;
  Parser *params;
  std::vector<spline> pathway;
  std::vector<double> pathway_r;
  std::string pafi_package;
  bool s_flag,has_pafi,cubic_spline;
private:
  /* nothing */
};


#endif
