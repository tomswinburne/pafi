#ifndef CSIM_H
#define CSIM_H

#include <mpi.h>

#include "LAMMPSSimulator.hpp"

class CustomSimulator : public LAMMPSSimulator {

public:
  CustomSimulator(MPI_Comm &instance_comm, Parser &p, int t) :
    LAMMPSSimulator(instance_comm, p, t){};

  void constrained_average(std::string SampleSteps) {
    std::string cmd;

    //if(local_rank==0) std::cout<<cmd<<std::endl;
    /*
      Here I can set up averages...
    */

    // run sampling
    cmd = "run "+SampleSteps;
    LAMMPSSimulator::run_commands(cmd);
    /*
      Here I can add them to results...
    */
  }

  int out_width = 16;


  // here I can print them
  void screen_output_header(double T) {
    LAMMPSSimulator::screen_output_header(T,out_width,false);
    std::cout<<std::setw(out_width)<<"Opt. \n";
  };

  void screen_output_line(double r){
    LAMMPSSimulator::screen_output_line(r,out_width,false);
    std::cout<<std::setw(out_width)<<"0.0 \n";
  };

  /*
    For the moment, we must "know" what is coming from LAMMPSSimulator
    results["MinEnergy"] = refE;
    results["preT"] = sampleT;
    results["postT"] = sampleT;
    results["aveF"] = *lmp_ptr * norm_mag;
    results["aveFstd"]=0.0; // from single worker
    results["stdF"] = *lmp_ptr * norm_mag * norm_mag;
    results["avePsi"] = *lmp_ptr;
    results["dXTangent"] = *lmp_ptr;
    results["MaxJump"] = sqrt(max_disp);
    results["Valid"] = double(bool(results["MaxJump"]<params->maxjump_thresh));
    results["MaxDev"] = dm;
  */
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

};

#endif
