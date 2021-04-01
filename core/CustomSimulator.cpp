#include "CustomSimulator.hpp"

/*
  Workflow:
  1) Load in configuration
  2) Apply hyperplane fixes
  3) Thermalize
  4) Declare PAFI free energy gradient and structural deviation time averages
  5) constrained_average() function (default : run SampleSteps)
  6) Extract time averages and add to results, compare structures

  LAMMPSSimulator results vector has :
    "MinEnergy", "preT", "postT", "aveF" (PAFI force), "stdF",
    "avePsi", "dXTangent", "MaxJump", "Valid", "MaxDev"
  */


CustomSimulator::CustomSimulator(MPI_Comm &instance_comm, Parser &p, int t)
  : LAMMPSSimulator(instance_comm, p, t) {
    /*
      anything extra to initialize ?
    */
    log_fields.push_back(std::make_pair("aveF",false));
  };

void CustomSimulator::sample(std::map<std::string,double> params, double *dev) {
  LAMMPSSimulator::run_commands("clear");
  made_fix=false; made_compute=false;
  LAMMPSSimulator::run_script("Input");
  LAMMPSSimulator::sample(params,dev);
};

//void CustomSimulator::populate() {


void CustomSimulator::constrained_average(std::string SampleSteps) {
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
