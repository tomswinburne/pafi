#ifndef LSIM_H
#define LSIM_H
#include <mpi.h>

#include "GeneralSimulator.hpp"

#include "lammps/lammps.h"
#include "lammps/library.h"

using namespace LAMMPS_NS;

class LAMMPSSimulator : public GeneralSimulator {

public:
  LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int t, int nr);

  /*
    Load xyz configuration from data file and put in vector
  */
  void load_config(std::string file_string,double *x);
  /*
    Parse and run script from configuration file
  */
  void run_script(std::string sn);
  /*
    Parse and run script from string with linebreaks
  */
  void run_commands(std::vector<std::string> strv);
  void run_commands(std::string strv);

  /*
    Fill configuration, path, tangent and tangent gradient. Return tangent norm
  */
  void populate(double r, double &norm_mag, double T);
  /*
    Rescale simulation cell
  */
  void rescale_cell(double T);
  /*
    Main sample run. Results vector should have thermalization temperature,
    sample temperature <f>, <f^2>, <psi> and <x-u>.n
  */
  void sample(double r, double T, std::map<std::string,double> &results,
                                                                  double *dev);

  std::string header(double mass);

  void lammps_write_dev(std::string fn, double r, double *dev, double *dev_sq);
  void lammps_dump_path(std::string fn, double r);

  double getEnergy();
  double getForceEnergy(double *f);

  // Fill 9D array with Lx, Ly, Lz, xy, xz, yz, then periodicity in x, y, z
  std::array<double,9> getCellData();

  void close();

  int error_count;

  std::string last_error();

  void log_error(std::string lc);
  void log_error(std::vector<std::string> lc);

  void gather(std::string name, int c, double *v);
  void gather(std::string name, int c, int *v);
  void scatter(std::string name, int c, double *v);
  void scatter(std::string name, int c, int *v);


  // LAMMPS specific
  int *species,*q, *image, *id;
private:
  int local_rank;
  LAMMPS *lmp;
  bool made_fix,made_compute;
  std::string last_error_message, last_command;

};

#endif
