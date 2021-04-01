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

// here I can print them
void CustomSimulator::screen_output_header(double T) {
  LAMMPSSimulator::screen_output_header(T,out_width,false);
  std::cout<<std::setw(out_width)<<"Opt. \n";
};

void CustomSimulator::screen_output_line(double r){
  LAMMPSSimulator::screen_output_line(r,out_width,false);
  std::cout<<std::setw(out_width)<<"0.0 \n";
};


/*
void fill_results(double r,double *ens_data) {
  LAMMPSSimulator::fill_results(r,ens_data,false);
  // this determines the splines. r is already added
  data_log.push_back(results["aveF"]); // splines[0]
  data_log.push_back(results["aveFstd"]); // splines[1]
  data_log.push_back(results["avePsi"]); // splines[2]
};

void integrate(std::string res_file, double &barrier){

  LAMMPSSimulator::integrate(res_file,barrier,false);

  std::ofstream out;
  out.open(res_file.c_str(),std::ofstream::out);

  out<<"# r av(F(r)) std(F(r)) ave(Psi)"<<std::endl;

  // choose integration parameters
  double diff_r = sample_r[sample_r.size()-1] - sample_r[0];
  double dr = diff_r / sample_r.size() / 10.0;
  double f_bar=0.0,ef_bar=0.0,f_bar_max=0.0;

  for(double r=sample_r[0];r<=sample_r[0]+diff_r;r+=dr) {

    f_bar -= dr/2.0 * splines[0](r);
    ef_bar += dr/2.0 * splines[1](r);
    f_bar_max = std::max(f_bar,f_bar_max);

    out<<r<<" "<<f_bar<<" "<<ef_bar<<" "<<splines[2](r)<<std::endl;

    f_bar -= dr/2.0 * splines[0](r);
    ef_bar += dr/2.0 * splines[1](r);
    f_bar_max = std::max(f_bar,f_bar_max);
  }

  barrier=f_bar_max;
};
*/
