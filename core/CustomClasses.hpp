#ifndef CSIM_H
#define CSIM_H

#include "GeneralGatherer.hpp"

#include "LAMMPSSimulator.hpp"

class CustomSimulator : public LAMMPSSimulator {

public:
  CustomSimulator(MPI_Comm &instance_comm, Parser &p, int t)
  : LAMMPSSimulator(instance_comm, p, t) {
    /*
      anything extra to initialize ?
    */
    log_fields.push_back(std::make_pair("aveF",false));
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

  void sample(std::map<std::string,double> params, double *dev) {
    LAMMPSSimulator::run_commands("clear");
    made_fix=false; made_compute=false;
    LAMMPSSimulator::run_script("Input");
    LAMMPSSimulator::sample(params,dev);
  };
};

class CustomGatherer : public GenericGatherer {
public:
  void screen_output_header() {
    GenericGatherer::screen_output_header(false);
    if(rank==0){
      std::cout<<std::setw(fw)<<" Opt";
      std::cout<<std::endl;
    }
  }
  void screen_output_line() {
    GenericGatherer::screen_output_line(false);
    if(rank==0){
      std::cout<<std::setw(fw)<<" 0.0";
      std::cout<<std::endl;
    }
  }
};

#endif
