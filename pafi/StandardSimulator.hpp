#ifndef CSIM_H
#define CSIM_H
#include "LAMMPSSimulator.hpp"

/*

*/

class StandardSimulator : public LAMMPSSimulator {

public:
  StandardSimulator(MPI_Comm &instance_comm, Parser &p, Holder &h, int t)
  : LAMMPSSimulator(instance_comm, p, h, t) {};

  void sample(Holder params, double *dev) {
    // reestablish
    LAMMPSSimulator::reset();
    // reload in case Input script has parameter dependence
    LAMMPSSimulator::run_script("Input");
    // Ensure we can access LAMMPS data
    LAMMPSSimulator::fill_lammps_vectors();
    // thermalize then execute constrained_average() function
    LAMMPSSimulator::sample(params,dev);
  };

  /*
    The main sample run fills a std::map<std::string,double> results with:
    0: MaxDev : total |X_0-X_t|
    1: MaxJump : max per-atom entry of |X_0-X_t|
    2: MinEnergy : potential energy following in-plane minimization
    3: Valid : considered valid so far: MaxJump< MaxJump in config file
    4: aveF : time averaged PAFI free energy gradient
    6: avePsi : time averaged PAFI MEP,MFEP overlap
    7: dXTangent : projection of [X_t-X_0] on plane normal. should be 0 ~ 1e-15
    8: postT : equipartition temperature following sample run
    9: preT : equipartition temperature before sample run
    10: stdF : Time variance of PAFI free energy gradient
  */
  void constrained_average() {
    // Just perform SampleSteps of constrained MD
    std::string cmd = "run "+parser->configuration["SampleSteps"];
    LAMMPSSimulator::run_commands(cmd);
  };

};


#endif
