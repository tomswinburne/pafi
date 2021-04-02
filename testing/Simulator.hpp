#ifndef CSIM_H
#define CSIM_H

#include "LAMMPSSimulator.hpp"

class LAMMPS_TI_Simulator : public LAMMPSSimulator {

public:
  LAMMPS_TI_Simulator(MPI_Comm &instance_comm, Parser &p, Holder &h, int t)
  : LAMMPSSimulator(instance_comm, p, h, t) {
    /*
      anything extra to initialize ?
    */
    std::cout<<"HHERERE"<<std::endl;
  };

  void constrained_average() {

    int nfix = lmp->modify->nfix;
    std::cout<<nfix<<" -> ";
    LAMMPSSimulator::run_script("PreSample");
    int pnfix = lmp->modify->nfix;
    std::cout<<pnfix<<std::endl;

    std::string cmd = "run "+SampleSteps;
    LAMMPSSimulator::run_commands(cmd);
    LAMMPSSimulator::run_script("PostSample");
  };

  void reset() {
    LAMMPSSimulator::reset();
  };

  void sample(Holder params, double *dev) {
    reset();

    // change potential....
    LAMMPSSimulator::run_script("Input");
    // set up natoms etc
    LAMMPSSimulator::fill_lammps_vectors();

    LAMMPSSimulator::run_script("Input");

    LAMMPSSimulator::sample(params,dev);
  };
};


#endif
