#ifndef CSIM_H
#define CSIM_H

#include "GeneralGatherer.hpp"

#include "LAMMPSSimulator.hpp"

class CustomSimulator : public LAMMPSSimulator {

public:
  CustomSimulator(MPI_Comm &instance_comm, Parser &p, int t);

  void constrained_average(std::string SampleSteps);

  void sample(std::map<std::string,double> params, double *dev);
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
};

#endif
