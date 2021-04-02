#ifndef CSIM_H
#define CSIM_H

#include "GeneralGatherer.hpp"

#include "LAMMPSSimulator.hpp"




class CustomGatherer : public GenericGatherer {
public:
  CustomGatherer(Parser &p, std::vector<double> pathway_r, int _nW, int di, int _rank) :
    GenericGatherer(p, pathway_r, _nW, di, _rank){};

  void screen_output_header(bool end=true) {
    GenericGatherer::screen_output_header(false);
    if(rank>0) return;
    std::cout<<std::setw(fw)<<" Opt"<<std::endl;
  };

  void screen_output_line(bool end=true) {
    GenericGatherer::screen_output_line(false);
    if(rank>0) return;
    std::cout<<std::setw(fw)<<" 0.0"<<std::endl;
  };

};


class CustomSimulator : public LAMMPSSimulator {

public:
  CustomSimulator(MPI_Comm &instance_comm, Parser &p, int t)
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
    LAMMPSSimulator::run_commands("clear");
    reset();
    LAMMPSSimulator::run_script("Input");
    LAMMPSSimulator::sample(params,dev);
  };
};


#endif
