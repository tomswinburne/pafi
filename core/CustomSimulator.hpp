#ifndef CSIM_H
#define CSIM_H


#include "LAMMPSSimulator.hpp"

class CustomSimulator : public LAMMPSSimulator {

public:
  CustomSimulator(MPI_Comm &instance_comm, Parser &p, int t);

  void constrained_average(std::string SampleSteps);

  void screen_output_header(double T);

  void screen_output_line(double r);
};

#endif
