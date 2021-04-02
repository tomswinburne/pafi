#ifndef CSIM_H
#define CSIM_H

#include "GenericGatherer.hpp"

#include "LAMMPSSimulator.hpp"




class TI_Gatherer : public GenericGatherer {
public:
  TI_Gatherer(Parser &p, std::vector<double> pathway_r, int _nW, int di, int _rank) :
    GenericGatherer(p, pathway_r, _nW, di, _rank){};

  // generated whenever a temperature cycle is complete
  void screen_output_header(bool end=true) {
    GenericGatherer::screen_output_header(false); // false : no endl;
    if(rank>0) return; // for flexibility
    std::cout<<std::setw(fw)<<" Opt"<<std::endl;
  };

  // generated whenever a reaction coordinate cycle is complete
  void screen_output_line(bool end=true) {
    GenericGatherer::screen_output_line(false); // false : no endl;
    if(rank>0) return; // for flexibility
    std::cout<<std::setw(fw)<<" 0.0"<<std::endl;
  };

};


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
    LAMMPSSimulator::run_script("Input");
    LAMMPSSimulator::sample(params,dev);
  };
};


#endif
