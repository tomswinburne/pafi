#ifndef LSIM_H
#define LSIM_H
#include <mpi.h>

#include "GeneralSimulator.hpp"

#include "lammps/lammps.h"
#include "lammps/library.h"
#include "lammps/atom.h"

class LAMMPSSimulator : public GeneralSimulator {

public:
  LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int rank);

  /*
    Load xyz configuration from data file and put in vector
  */
  void load_config(std::string file_string,std::vector<double> &x);

  /*
    Parse and run script from configuration file
  */
  void run_script(std::string sn);

  /*
    Parse and run script from string with linebreaks
  */
  void run_commands(std::string strv);

  /*
    Parse and run script from vector of strings
  */
  void run_commands(std::vector<std::string> strv);

  /*
    Fill configuration, path, tangent and tangent gradient. Return tangent norm
  */
  void populate(double r, double *scale, double &nm, bool justpos=false);

  /*
    Fill vec from locally stored knots
  */
  void fill_path(double r,int der,std::vector<double> &vec, double *scale);

  /*
    Rescale simulation cell with scale vector
  */
  void rescale_cell(double *scale);

  /*
    Main sample run. Results vector should have thermalization temperature,
    sample temperature <f>, <f^2>, <psi> and <x-u>.n
  */
  void sample(double r, double T, double *results, double *dev);

  /*
    Dump deviations
  */
  void write_dev(std::string fn, double *dev, double *dev_sq);

  /*
    Return Energy
  */
  double getEnergy();

  /*
    Return Force and Energy
  */
  double getForceEnergy(double *f);


  // Fill 9D array with Lx, Ly, Lz, xy, xz, yz, then periodicity in x, y, z
  std::array<double,9> getCellData();

  void close();

  double temperature;
private:
  LAMMPS_NS::LAMMPS *lmp=NULL;
  //std::vector<int> id, species, image;
  //std::vector<double> q;
};

#endif
