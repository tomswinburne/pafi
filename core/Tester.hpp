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
    for(i=0;i<g.dsize;i++) g.all_data[i]=g.data[i]; // hack for test
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
    double dr = diff_r / sample_r.size() / 10.0;
    double F_bar = 0., E_bar=0., f=0.;
    for(auto e: dE) E_bar = std::max(E_bar,e-dE[0]);
    for(double r=sample_r[0];r<=sample_r[0]+diff_r;r+=dr) {
      f -= dr/2.0 * Fspl(r);
      F_bar = std::max(F_bar,f);
      f -= dr/2.0 * Fspl(r);
    }
    std::cout<<"\n\n******************************************************************\n\n";
    std::cout<<"\tEnergy Barrier ~= "<<E_bar<<"eV, Force Integration Barrier ~= "<<F_bar<<" eV\n";
    std::cout<<"\tONLY SIMPLE TESTS PERFORMED HERE!\n"
    "\tPLEASE USE pafi-path-test TO TEST PATHWAY FOR FORCE INTEGRATION\n"
    "\n******************************************************************\n\n"<<std::endl;

  }

  // close down LAMMPS instances
  sim.close();

  MPI_Barrier(world);

  // close down MPI
  MPI_Comm_free(&instance_comm);
};
