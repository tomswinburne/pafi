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
};

#endif
