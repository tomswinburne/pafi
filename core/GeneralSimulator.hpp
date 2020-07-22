#ifndef SIM_H
#define SIM_H

#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <mpi.h>

#include "Constants.hpp"
#include "Parser.hpp"
#include "Boundary.hpp"
#include "Spline.hpp"

class GeneralSimulator {

public:
  GeneralSimulator(MPI_Comm &instance_comm, Parser &p, int rank);
  
  void expansion(double T,double *scale);

  void make_path(std::vector<std::string> knot_list);

  void fill_path(double r,int der,std::vector<double> &vec);

  int getNatoms(){
    return natoms;
  };

  // LAMMPS DEPENDENT

  virtual void load_config(std::string file_string,std::vector<double> &x){};

  virtual void run_script(std::string sn){};

  virtual void run_commands(std::vector<std::string> strv){};

  virtual void setup(double r, double T){};

  virtual void sample(double r, double T, double *results, double *dev){};

  virtual double getEnergy(){};

  virtual double getForceEnergy(double *f){};

  virtual std::array<double,9> getCellData(){};

  virtual void close(){};

  virtual void populate(double r){};

  virtual void rescale_cell(double *scale){};

  virtual void sample(double r, double T, std::vector<double> &results){};

  virtual void write_dev(std::string fn, double r, double *dev, double *dev_sq);

protected:
  double refE, position;
  double norm_com[3], norm_mag;
  int natoms, tag, nknots, nlocal, local_rank, local_size;
  MinImage pbc;
  Parser *params;
  std::vector<spline> pathway;
  std::vector<int> loc_id;
  std::vector<double> loc_data, glo_data, t;
  MPI_Comm *comm;
private:
  /* nothing */
};


#endif
