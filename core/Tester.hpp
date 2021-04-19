#include <iostream>
#include <sstream>
#include <iomanip>

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <array>
#include <map>
#include <list>
#include <cmath>

// PAFI files
#include "ConstantsTypes.hpp"
#include "Parser.hpp"
#include "Spline.hpp"

template <class SimulatorTemplate,class GathererTemplate>
void test(MPI_Comm &world,std::string parser_file,bool lammps_prep) {

  int rank, nProcs, i, j;

  MPI_Comm_rank(world,&rank);
  MPI_Comm_size(world,&nProcs);

  // ************************ READ CONFIG FILE **********************************
  Parser parser(parser_file,false);

  if(!parser.xml_success) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Configuration file could not be read!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }
  parser.overwrite_xml(nProcs);
  parser.set_parameters();
  // ************************ READ CONFIG FILE ***********************************


  // ************************ DUMP FOLDER *********************************
  int dump_index=-1;
  if(rank==0) parser.find_dump_file(dump_index);
  MPI_Bcast(&dump_index,1,MPI_INT,0,world);
  MPI_Barrier(world);
  if(dump_index<0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"PAFI could not write to output path : does the directory ";
      std::cout<<parser.dump_dir<<" exist? Exiting!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }
  // ************************ DUMP FOLDER *********************************



  // ******************* SET UP WORKERS ***************************************
  const int nWorkers = 1;
  const int instance = 0;
  const int local_rank = rank;
  const int min_valid = 0;

  // LAMMPS communicators
  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  // see GlobalSeed
  parser.seed(instance);

  // set up data gatherer
  GathererTemplate g(parser,nWorkers,dump_index,rank);
  MPI_Bcast(&(g.initialized),1,MPI_INT,0,world);
  if(g.initialized==0) exit(-1);

  SimulatorTemplate sim(instance_comm,parser,g.params,instance);
  if(!sim.has_pafi) exit(-1);
  if(rank==0)  std::cout<<parser.welcome_message();

  sim.make_path(parser.PathwayConfigurations);
  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";
  g.special_r_overwrite(sim.pathway_r);

  // ******************* SET UP WORKERS ****************************************


  // ********************** SAMPLING *******************************************

  // generic - deviation
  const int vsize = 3 * sim.natoms;
  double *dev = new double[vsize*(1+(rank==0))];

  int valid[2] = {1,1};

  if(rank==0) {
    std::cout<<"\n\n\n*****************************\n\n\n";
    std::cout<<"PAFI TEST RUN \n";
    std::cout<<"\n\n\n*****************************\n\n\n";
  }


  g.screen_output_header();

  int fileindex = 0;
  std::string dump_file_name;

  while(true) {
    // sample
    for(i=0;i<vsize;i++) dev[i] = 0.0;
    sim.sample(g.params, dev); // sim*(parser, dev)

    g.prepare(sim.results); // allocate memory if not already done
    if(rank==0) for(i=0;i<g.dsize;i++) g.all_data[i]=g.data[i]; // hack for test
    int total_valid = g.collate(valid+1);

    if(lammps_prep) {
      dump_file_name = "dumps/pafipath."+std::to_string(fileindex)+".dat";
      sim.lammps_dump_path(dump_file_name,g.params["ReactionCoordinate"]);
      fileindex++;
    } else {
      if(rank==0 and parser.write_dev) {
        sim.write_dev(g.dev_file,g.params["ReactionCoordinate"],dev+vsize);
      }
    }

    g.next(); // wipe ens_data
    MPI_Barrier(world);

    if(g.finished()) break;
  }
  if(rank==0) {
    g.close(); // just the dump files
    std::vector<double> dF,maxjumpr,dE,sample_r;
    for(auto ens_res : g.all_ens_results) {
      Holder p = ens_res.first;
      auto res = ens_res.second;
      sample_r.push_back(p["ReactionCoordinate"]);
      dF.push_back(res["aveF"].first);
      dE.push_back(res["MinEnergy"].first);
      maxjumpr.push_back(res["MaxJump"].first);
    }
    spline Fspl;
    Fspl.set_points(sample_r,dF);
    double diff_r = sample_r[sample_r.size()-1] - sample_r[0];
    double dr = diff_r / sample_r.size() / 20.0;
    double F_bar = 0., E_bar=0., f=0.;

    for(auto e: dE) E_bar = std::max(E_bar,e-dE[0]);
    for(double r=sample_r[0];r<=sample_r[0]+diff_r;r+=dr) {
      f += dr/2.0 * Fspl(r);
      F_bar = std::max(F_bar,f);
      f += dr/2.0 * Fspl(r);
    }
    std::cout<<"\n\n******************************************************************\n\n";

    if(lammps_prep) {
      std::cout<<"\tEnergy Barrier ~= "<<E_bar<<"eV, Force integration Barrier ~= "<<F_bar<<" eV\n";
      std::cout<<"\tONLY SIMPLE TESTS PERFORMED HERE!\n"
      "\tPLEASE USE pafi-path-test TO TEST PATHWAY FOR FORCE INTEGRATION\n"
      "\n******************************************************************\n\n"<<std::endl;
    } else {

      std::cout<<"Absolute Forces and differences between knots: \n";


      double _dF,_ddF,ddF[2] = {0.,0.},F=0.;
      bool warning=false;
      int i=0;
      for (;i<dF.size()-1;i++) {
        _ddF = std::fabs(dF[i+1]-dF[i]);
        _dF = -(dF[i]+dF[i+1])/2. * (sample_r[i+1]-sample_r[i]);

        std::cout<<"\n\t Knot "<<i+1<<": r = "<<sample_r[i]<<" |dF| = "<<std::fabs(dF[i])<<" , F = "<<F<<"\n";
        std::cout<<"\n\t\t   Max | X_postmin - X | = "<<maxjumpr[i];
        if(maxjumpr[i]>0.02) {
          warning = true;
          std::cout<<" !! This should be zero for a properly discretized minimum energy path!\n"
          " Values greater than ~0.02 should be considered risky;"
          " Perhaps consider more integration points here. WARNING";
        }

        std::cout<<"\n\t\t        | dF_next - dF | = "<<_ddF;
        std::cout<<"\n\t\t Approx |  F_next - F  | = "<<_dF;
        std::cout<<"\n\t\t        |  r_next - r  | = "<<sample_r[i+1]-sample_r[i]<<"\n";

        if(_ddF<0.02 && std::fabs(dF[i])<0.02) {
          warning = true;
          std::cout<<"Very flat segment! "
          "Perhaps remove knot "<<i+1<<" or "<<i+2<<"?\n. FAIL";
        }
        double rel_dF = (_dF/F_bar*sample_r.size()/2.0); // =1 for triangle barrier
        if(rel_dF>5.0) {
          warning = true;
          std::cout<<"Large free energy change! "
          "Consider more integration points here. FAIL";
        }

        ddF[0] += _ddF/double(dF.size()-1);
        ddF[1] += _ddF*_ddF/double(dF.size()-1);

        F += -(dF[i]+dF[i+1])/2. * (sample_r[i+1]-sample_r[i]);
        std::cout<<"\n -------- \n";
      }
      ddF[1] -= ddF[0]*ddF[0];
      std::cout<<"\n\t Knot "<<i+1<<": r = "<<sample_r[i]<<" |dF| = "<<std::fabs(dF[i])<<" , F = "<<F<<"\n";
      std::cout<<"\n\t\t   Max | X_postmin - X | = "<<maxjumpr[i];
      std::cout<<"\n -------- \n";
      std::cout<<"\n\n\tAverage, Std in d|dF|:"<<ddF[0]<<" , "<<sqrt(ddF[1])<<std::endl;



      std::cout<<"\n\n--------------- integration checks at zero temperature -----\n";
      std::cout<<"\n\tEnergy Barrier: "<<std::setprecision(5)<<E_bar;
      std::cout<<" eV, Force integration Barrier: "<<std::setprecision(5)<<F_bar<<" eV";
      std::cout<<"\n\n\tAbsolute error of: "<<std::setprecision(5)<<std::fabs(F_bar-E_bar)*1000.<<" meV ";
      double pc_prec = (F_bar/E_bar-1.0)*100.;
      std::cout<<" ("<<std::setprecision(5)<<pc_prec<<"%)"<<std::endl;
      if(std::fabs(pc_prec)<=2) {
        std::cout<<"\n\tError < 2% - within likely potential accuracy, OK for sampling!\n";
      } else {
        warning=true;
        std::cout<<"\n\tError > 2%, could be too high...\n";
      }
      if(std::fabs(pc_prec)>=0.5) {
        std::cout<<"\t\n\tTo reduce this error: ";
        std::cout<<"\n\t\t More NEB images and/or increasing nPlanes in config.xml\n\n"<<std::endl;
      }

      if(warning) std::cout<<"\nPathway checks failed / warnings generated! See above\n"<<std::endl;
      else std::cout<<"\nPathway checks passed!\n"<<std::endl;
    }


  }

  // close down LAMMPS instances
  sim.close();

  MPI_Barrier(world);

  // close down MPI
  MPI_Comm_free(&instance_comm);
};
