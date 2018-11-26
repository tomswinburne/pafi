#ifndef LSIM_H
#define LSIM_H
#include <mpi.h>

#include "GeneralSimulator.hpp"

#include "lammps/lammps.h"
#include "lammps/library.h"

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
  void run_commands(std::vector<std::string> strv);
  /*
    Fill configuration, path, tangent and tangent gradient. Return tangent norm
  */
  void populate(double r, double scale, double &norm_mag);
  /*
    Rescale simulation cell
  */
  void rescale_cell(double scale);
  /*
    Main sample run. Results vector should have thermalization temperature,
    sample temperature <f>, <f^2>, <psi> and <x-u>.n
  */
  void sample(double r, double T, double *results, double *dev);

  std::string header();

  void lammps_dev_write(std::string fn, double r, double *dev, double *dev_sq);

  void lammps_path_write(std::string fn, double r);

  double getEnergy();

  // Fill 9D array with Lx, Ly, Lz, xy, xz, yz, then periodicity in x, y, z
  std::array<double,9> getCellData();

  void close();

  // LAMMPS specific
  std::vector<int> species,q,image,id;
private:
  void *lmp;
};

#endif
