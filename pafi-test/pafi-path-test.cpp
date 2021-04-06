#include "pafi-test.hpp"

int main(int narg, char **arg) {

  MPI_Init(&narg,&arg);
  int rank, nProcs, i, j;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // ************************ READ CONFIG FILE **********************************
  Parser parser("./config.xml",true);

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
  MPI_Bcast(&dump_index,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(dump_index<0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Could not write to output path / find directory! Exiting!\n";
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
  Gatherer g(parser,nWorkers,dump_index,rank);
  MPI_Bcast(&(g.initialized),1,MPI_INT,0,MPI_COMM_WORLD);
  if(g.initialized==0) exit(-1);



  Simulator sim(instance_comm,parser,g.params,instance);
  if(!sim.has_pafi) exit(-1);
  if(rank==0)  std::cout<<parser.welcome_message();

  sim.make_path(parser.PathwayConfigurations);
  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";
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

  while(true) {
    // sample
    for(i=0;i<vsize;i++) dev[i] = 0.0;
    sim.sample(g.params, dev); // sim*(parser, dev)

    g.prepare(sim.results); // allocate memory if not already done
    for(i=0;i<g.dsize;i++) g.all_data[i]=g.data[i]; // hack for test
    int total_valid = g.collate(valid+1);

    if(rank==0 and parser.write_dev)
      sim.write_dev(g.dev_file,g.params["ReactionCoordinate"],dev+vsize);

    g.next(); // wipe ens_data
    MPI_Barrier(MPI_COMM_WORLD);

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

    std::cout<<"\nT=0K test run complete, testing path....\n\n";
    std::cout<<"Absolute Forces and differences between knots: \n";

    double _dF,_ddF,ddF[2] = {0.,0.},_F=0.;
    bool warning=false;
    for (i=0;i<dF.size()-1;i++) {
      _ddF = std::fabs(dF[i+1]-dF[i]);
      _dF = -(dF[i]+dF[i+1])/2. * (sample_r[i+1]-sample_r[i]);
      std::cout<<"\n\t Knot "<<i+1<<": r = "<<sample_r[i]<<" |dF| = "<<std::fabs(dF[i])<<" , F = "<<_F<<", E = "<<dE[i]-dE[0]<<"\n";
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

      if(_dF>0.1) {
        warning = true;
        std::cout<<"Large free energy change! "
        "Consider more integration points here. FAIL";
      }

      ddF[0] += _ddF/double(dF.size()-1);
      ddF[1] += _ddF*_ddF/double(dF.size()-1);

      _F += -(dF[i]+dF[i+1])/2. * (sample_r[i+1]-sample_r[i]);
      std::cout<<"\n -------- \n";
    }
    ddF[1] -= ddF[0]*ddF[0];
    std::cout<<"\n\t Knot "<<i+1<<": r = "<<sample_r[i]<<" |dF| = "<<std::fabs(dF[i])<<" , F = "<<_F<<"\n";
    std::cout<<"\n\t\t   Max | X_postmin - X | = "<<maxjumpr[i];
    std::cout<<"\n -------- \n";
    std::cout<<"\n\n\tAverage, Std in d|dF|:"<<ddF[0]<<" , "<<sqrt(ddF[1])<<std::endl;



    std::cout<<"\n\n--------------- Integration checks at zero temperature -----\n";
    std::cout<<"\n\tEnergy Barrier: "<<std::setprecision(5)<<E_bar;
    std::cout<<" eV, Force Integration Barrier: "<<std::setprecision(5)<<F_bar<<" eV";
    std::cout<<"\n\n\tAbsolute error of: "<<std::setprecision(5)<<std::fabs(F_bar-E_bar)*1000.;
    std::cout<<" meV "<<" ("<<std::setprecision(5)<<(F_bar/E_bar-1.0)*100.<<"%)";
    if(std::fabs(F_bar-E_bar)<=0.002) {
      std::cout<<"  ** Even if small, this error should be added to PAFI error bars **"<<std::endl;
      std::cout<<"\n\tError < 2meV - within likely potential accuracy, great for sampling! OK\n";
    } else if(std::fabs(F_bar-E_bar)<0.01) {
      std::cout<<"  ** This error should be added to PAFI error bars **"<<std::endl;
      std::cout<<"\n\tAbove target error of ~2meV, but < 10meV probably OK for sampling CAUTION\n";
    } else {
      warning=true;
      std::cout<<"\n\tError > 10meV  - if using this setup, the "
      "error of "<<F_bar-E_bar<<" eV MUST be added to"
      "PAFI error bars, and could be significant! FAIL\n";
    }
    if(std::fabs(F_bar-E_bar)>=0.001) {
      std::cout<<"\n\tTo reduce this error: ";
      std::cout<<"\n\t More NEB images and/or increasing nPlanes in config.xml";
      std::cout<<"\n\t Try Rediscretize==0/1 to change how images are interpolated\n\n"<<std::endl;
    }

    if(warning) std::cout<<"\nPathway checks failed / warnings generated! See above\n"<<std::endl;
    else std::cout<<"\nPathway checks passed!\n"<<std::endl;
  }

  // close down LAMMPS instances
  sim.close();

  MPI_Barrier(MPI_COMM_WORLD);

  // close down MPI
  MPI_Comm_free(&instance_comm);


  MPI_Finalize();

}
