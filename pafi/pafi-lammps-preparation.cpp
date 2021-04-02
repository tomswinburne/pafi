#include "pafi.hpp"

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

  Simulator sim(instance_comm,parser,instance);
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

  // set up data gatherer
  Gatherer g(parser,sim.pathway_r,nWorkers,dump_index,rank);
  MPI_Bcast(&(g.initialized),1,MPI_INT,0,MPI_COMM_WORLD);
  if(g.initialized==0) exit(-1);

  g.screen_output_header();
  int fileindex = 0;
  while(true) {
    // sample
    for(i=0;i<vsize;i++) dev[i] = 0.0;
    sim.sample(g.params, dev); // sim*(parser, dev)

    g.prepare(sim.results); // allocate memory if not already done
    for(i=0;i<g.dsize;i++) g.all_data[i]=g.data[i]; // hack for test
    int total_valid = g.collate(valid+1);

    double r = g.params["ReactionCoordinate"];
    std::string dump_file_name = "dumps/pafipath."+std::to_string(fileindex)+".dat";
    sim.lammps_dump_path(dump_file_name,r);

    fileindex++;

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
    std::cout<<"\n\n******************************************************************\n\n";
    std::cout<<"\tEnergy Barrier ~= "<<E_bar<<"eV, Force Integration Barrier ~= "<<F_bar<<" eV\n";
    std::cout<<"\tONLY SIMPLE TESTS PERFORMED HERE!\n"
    "\tPLEASE USE pafi-path-test TO TEST PATHWAY FOR FORCE INTEGRATION\n"
    "\n******************************************************************\n\n"<<std::endl;

  }

  // close down LAMMPS instances
  sim.close();

  MPI_Barrier(MPI_COMM_WORLD);

  // close down MPI
  MPI_Comm_free(&instance_comm);


  MPI_Finalize();

}
