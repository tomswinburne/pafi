#ifndef LSIM_H
#define LSIM_H

#include "GeneralSimulator.hpp"

#include "lammps/lammps.h"
#include "lammps/library.h"

using namespace LAMMPS_NS;

class LAMMPSSimulator : public GeneralSimulator {

public:
  LAMMPSSimulator(MPI_Comm &instance_comm, Parser &p, int t);

  /*
    Load xyz configuration from data file and put in vector
  */
  virtual void load_config(std::string file_string,double *x);
  /*
    Parse and run script from configuration file
  */
  virtual void run_script(std::string sn);
  /*
    Parse and run script from string with linebreaks
  */
  virtual void run_commands(std::vector<std::string> strv);
  virtual void run_commands(std::string strv);

  /*
    Fill configuration, path, tangent and tangent gradient. Return tangent norm
  */
  virtual void populate(double r, double &norm_mag, double T);
  /*
    Rescale simulation cell
  */
  virtual void rescale_cell(double T);
  /*
    Main sample run. Results vector should have thermalization temperature,
    sample temperature <f>, <f^2>, <psi> and <x-u>.n
  */
  virtual void sample(double r, double T, double *dev);

  virtual void constrained_average(std::string SampleSteps);

  virtual void screen_output_header(double T, int fw=18, bool end=true);
  virtual void screen_output_line(double r, int fw=18, bool end=true);

  virtual void fill_results(double r,double *ens_data, bool end=true);
  virtual void end_of_cycle(std::string res_file,bool end=true);

  virtual std::string header(double mass);

  virtual void lammps_write_dev(std::string fn, double r, double *dev, double *dev_sq);
  virtual void lammps_dump_path(std::string fn, double r);

  virtual double getEnergy();
  virtual double getForceEnergy(double *f);
  virtual double get_fix(std::string fixid,int type,int index);


  virtual void close();

  virtual std::string last_error();

  virtual void log_error(std::string lc);
  virtual void log_error(std::vector<std::string> lc);

  virtual void gather(std::string name, int c, double *v);
  virtual void gather(std::string name, int c, int *v);
  virtual void scatter(std::string name, int c, double *v);
  virtual void scatter(std::string name, int c, int *v);

  // Fill 9D array with Lx, Ly, Lz, xy, xz, yz, then periodicity in x, y, z
  std::array<double,9> getCellData();
  // LAMMPS specific
  int *species,*q, *image, *id;
  double *lt; // for scattering
protected:
  LAMMPS *lmp;
  bool made_fix,made_compute;
  std::string last_command;

};

#endif
