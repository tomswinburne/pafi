#ifndef CSIM_H
#define CSIM_H

#include "LAMMPSSimulator.hpp"

class LAMMPS_TI_Simulator : public LAMMPSSimulator {

public:
  LAMMPS_TI_Simulator(MPI_Comm &instance_comm, Parser &p, int t)
  : LAMMPSSimulator(instance_comm, p, t) {
    /*
      anything extra to initialize ?
    */
  };

  void constrained_average(std::string SampleSteps) {
    std::string cmd;
    /*
      Set up fixes + time averages...
    */
    // run sampling
    cmd = "run "+SampleSteps;
    LAMMPSSimulator::run_commands(cmd);
    /*
      Extract and add to results...
    */
  };

  void reset() {
    LAMMPSSimulator::reset();
  };

  void sample(Holder params, double *dev) {
    reset();
    auto commands = parser->Script("input");
    // change potential....
    LAMMPSSimulator::run_script("Input");
    // set up natoms etc
    LAMMPSSimulator::fill_lammps_vectors();

    LAMMPSSimulator::run_script("Input");

    LAMMPSSimulator::sample(params,dev);
  };
};


#endif
