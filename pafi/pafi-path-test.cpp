#include "pafi.hpp"

int main(int narg, char **arg) {

  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Load input file
  Parser params("./config.xml",true);

  if(!params.xml_success) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Configuration file could not be read!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }

  params.overwrite_xml(nProcs);
  params.set_parameters();

  // check for dump folder
  std::ofstream raw;
  int dump_index=-1;
  if(rank==0) params.find_dump_file(raw,dump_index);
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

  if(rank==0)  {
    std::cout<<params.welcome_message();
    std::cout<<"\n\n\n*****************************\n\n\n";
    std::cout<<"TEST RUN AT T=0K FOR SPLINE CHECK\n";
    std::cout<<"\n\n\n*****************************\n\n\n";
  }

  MPI_Barrier(MPI_COMM_WORLD);

  const int nWorkers = 1;
  const int instance = 0;
  const int local_rank = rank;
  const int nRepeats = 1;
  const int nRes = 8; // TODO remove this

  params.seed(0);
  params.seed(123);
  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,0,0,&instance_comm);

  Simulator sim(instance_comm,params,0);
  if(!sim.has_pafi) exit(-1);

  sim.make_path(params.PathwayConfigurations);

  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";


  // Declare all variables

  const int vsize = 3 * sim.natoms;
  const int rsize = nRes*nWorkers;

  int *valid = new int[nWorkers];
  double *local_res = new double[rsize];
  double *local_dev = new double[vsize];
  double *local_dev_sq = new double[vsize];
  double *all_dev = NULL;
  double *all_dev_sq = NULL;
  double *all_res = NULL;

  if (rank == 0) {
    all_dev = new double[vsize];
    all_dev_sq = new double[vsize];
    all_res = new double[rsize];
  }

  int total_valid_data, totalRepeats, total_invalid_data;
  double temp, dr, t_max_jump, p_jump, E_init, E_max=0.0;
  std::string rstr,Tstr,dump_fn,fn,dump_suffix;
  std::map<std::string,double> results;

  std::vector<double> integr, dfer, dfere, psir, maxjumpr;
  std::vector<double> valid_res, invalid_res;
  // TODO : replace
  std::list<int> raw_dump_indicies = {2,3,4,7};// <f> std(f) <Psi> <Jump>



  if (params.nPlanes>1)
    dr = (params.stopr-params.startr)/(double)(params.nPlanes-1);
  else dr = 0.1;

  double T = 0.0;

  Tstr = std::to_string((int)(T));

  dump_suffix = "_"+Tstr+"K_"+std::to_string(dump_index);

  if(rank==0) {
    fn = params.dump_dir + "/raw_ensemble_output"+dump_suffix;
    std::cout<<"opening dump file "<<fn<<std::endl;
    raw.open(fn.c_str(),std::ofstream::out);
    std::cout<<"\nStarting T="+Tstr+"K run\n\n";
  }

  for(int i=0;i<nWorkers;i++) valid[i]=0;
  for(int i=0;i<rsize;i++) local_res[i]=0.0;
  for(int i=0;i<vsize;i++) local_dev[i] = 0.0;
  for(int i=0;i<vsize;i++) local_dev_sq[i] = 0.0;
  results.clear();
  integr.clear();
  dfer.clear();
  dfere.clear();
  psir.clear();
  valid_res.clear();
  invalid_res.clear();
  maxjumpr.clear();

  if (rank == 0) {
    std::cout<<"<> == time averages,  av/err over ensemble"<<std::endl;
    std::cout<<std::setw(5)<<"r";
    std::cout<<std::setw(20)<<"av(<Tpre>)";
    std::cout<<std::setw(20)<<"av(<Tpost>)";
    std::cout<<std::setw(20)<<"av(<dF/dr>)";
    std::cout<<std::setw(20)<<"err(<dF/dr>)";
    std::cout<<std::setw(20)<<"av(|(<X>-U).N|)";
    std::cout<<std::setw(20)<<"av(<N_true>.N)";
    std::cout<<std::setw(20)<<"Max Jump";
    std::cout<<std::setw(20)<<"P(Valid)";
    std::cout<<"\n";
  }

  for(auto r: sim.sample_r) {
    valid_res.clear();
    invalid_res.clear();
    rstr = std::to_string(r);
    for(int i=0;i<rsize;i++) local_res[i] = 0.0;
    for(int i=0;i<vsize;i++) local_dev[i] = 0.0;

    // store running average of local_dev in local_dev_sq
    total_valid_data=0;
    totalRepeats=0;
    t_max_jump=0.0;
    while (total_valid_data<=int(params.redo_thresh*nWorkers*nRepeats)) {

      sim.sample(r, T, results, local_dev);


      if(r==sim.sample_r[0]) E_init = results["MinEnergy"];

      E_max = std::max(results["MinEnergy"]-E_init,E_max);
      for(int i=0;i<vsize;i++) local_dev_sq[i] = local_dev[i]*local_dev[i];

      totalRepeats++;

      if(local_rank == 0) {
        if(sim.error_count>0) std::cout<<sim.last_error()<<std::endl;
        local_res[instance*nRes + 0] = results["preT"];
        local_res[instance*nRes + 1] = results["postT"];
        local_res[instance*nRes + 2] = results["aveF"];
        local_res[instance*nRes + 3] = results["stdF"];
        local_res[instance*nRes + 4] = results["aveP"];
        local_res[instance*nRes + 5] = results["TdX"];
        local_res[instance*nRes + 6] = results["MaxDev"];
        local_res[instance*nRes + 7] = results["MaxJump"];
      }

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(local_res,all_res,rsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      if (rank==0) {
        for(int i=0;i<nWorkers;i++) {
          t_max_jump = std::max(t_max_jump,all_res[i*nRes+7]);
          valid[i] = int(all_res[i*nRes+7]<params.maxjump_thresh);
          if(valid[i]) {
            for(int j=0;j<nRes;j++) valid_res.push_back(all_res[i*nRes+j]);
          } else {
            for(int j=0;j<nRes;j++) invalid_res.push_back(all_res[i*nRes+j]);
          }
        }
      }

      // Broadcast validity
      MPI_Bcast(valid,nWorkers,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);

      // nullify invalid batches
      if(valid[instance]==0 && !params.postDump) for(int i=0;i<vsize;i++) {
        local_dev[i]=0.0;
        local_dev_sq[i]=0.0;
      }

      // add all to total
      MPI_Reduce(local_dev,all_dev,vsize,
          MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(local_dev_sq,all_dev_sq,vsize,
          MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);

      for(int i=0;i<nWorkers;i++) total_valid_data += valid[i];
      p_jump = double(total_valid_data) / double(nWorkers*totalRepeats);

      if(totalRepeats>=params.maxExtraRepeats+nRepeats) break;

      if(p_jump < 0.1 ) {
        if(rank==0)
          std::cout<<"Reference path too unstable for sampling!"<<std::endl;
        break;
      }
    }
    // collate
    if(rank==0) {
      total_invalid_data = totalRepeats*nWorkers - total_valid_data;

      // deviation vectors
      if(params.postDump) for(int i=0;i<vsize;i++) {
        all_dev[i]/=double(totalRepeats*nWorkers);
        all_dev_sq[i]/=double(totalRepeats*nWorkers);
      } else if(total_valid_data>0) for(int i=0;i<vsize;i++) {
        all_dev[i]/=double(total_valid_data);
        all_dev_sq[i]/=double(total_valid_data);
      }

      //dump_fn = params.dump_dir+"/dev_"+rstr+dump_suffix+".dat";
      //sim.write_dev(dump_fn,r,all_dev,all_dev_sq);

      // raw output
      raw<<total_valid_data<<" ";
      for(auto j: raw_dump_indicies)
        for(int i=0;i<total_valid_data;i++) raw<<valid_res[i*nRes+j]<<" ";
      for(auto j: raw_dump_indicies)
        for(int i=0;i<total_invalid_data;i++) raw<<invalid_res[i*nRes+j]<<" ";
      raw<<std::endl;

      double *final_res = new double[2*nRes];
      for(int j=0;j<2*nRes;j++) final_res[j] = 0.;

      // average
      for(int i=0;i<total_valid_data;i++) for(int j=0;j<nRes;j++)
        final_res[j] += valid_res[i*nRes+j] / double(total_valid_data);

      for(int i=0;i<total_valid_data;i++) for(int j=0;j<nRes;j++) {
        temp = valid_res[i*nRes+j]-final_res[j];
        final_res[j+nRes] += temp * temp;
      }
      // under assumption that sample time is longer than autocorrelations,
      // expected variance in time average is ensemble variance / nWorkers
      // raw_output gives data to confirm this assumption (CLT with grouping)
      if(total_valid_data>0) for(int j=0;j<nRes;j++)
        final_res[j+nRes] = sqrt(final_res[j+nRes])/double(total_valid_data);

      integr.push_back(r);
      dfer.push_back(final_res[2]); // <dF/dr>
      dfere.push_back(final_res[2+nRes]);//"err(<dF/dr>)"
      psir.push_back(final_res[4]); // <Psi>
      maxjumpr.push_back(t_max_jump);
      std::cout<<std::setw(5)<<r;//"r"
      std::cout<<std::setw(20)<<final_res[0];//"av(<Tpre>)"
      std::cout<<std::setw(20)<<final_res[1];//"av(<Tpost>)"
      std::cout<<std::setw(20)<<final_res[2];//"av(<dF/dr>)"
      std::cout<<std::setw(20)<<final_res[2+nRes];//"err(<dF/dr>)"
      std::cout<<std::setw(20)<<final_res[5];//"av(|<X>-U).(dU/dr)|)"
      std::cout<<std::setw(20)<<final_res[4];//"av(Psi)"
      std::cout<<std::setw(20)<<t_max_jump;// max jump
      std::cout<<std::setw(20)<<p_jump;// ratio of jumps
      std::cout<<"\n";
    }
  }

  // free_energy_profile
  if(rank==0) {
    raw.close();
    std::cout<<"\nT="+Tstr+"K run complete, testing path....\n\n";

    std::cout<<"Absolute Forces and differences between knots: \n";


    spline dFspl;

    dFspl.set_points(integr,dfer);
    double dF_dense_int = 0.0,F_bar_dense=0.0,E_bar_dense=E_max;
    double int_dr = (sim.sample_r[sim.sample_r.size()-1]-sim.sample_r[0])/30.;

    for (double r = sim.sample_r[0]; r <= sim.sample_r[sim.sample_r.size()-1]; r += int_dr) {
      dF_dense_int -= dFspl(r)/2.0 * int_dr;
      F_bar_dense = std::max(F_bar_dense,dF_dense_int);
      dF_dense_int -= dFspl(r)/2.0 * int_dr;
    }


    double _dF,_ddF,ddF[2] = {0.,0.},F=0.;
    bool warning=false;
    int i=0;
    for (;i<dfer.size()-1;i++) {
      _ddF = std::fabs(dfer[i+1]-dfer[i]);
      _dF = -(dfer[i]+dfer[i+1])/2. * (integr[i+1]-integr[i]);

      std::cout<<"\n\t Knot "<<i+1<<": r = "<<integr[i]<<" |dF| = "<<std::fabs(dfer[i])<<" , F = "<<F<<"\n";
      std::cout<<"\n\t\t   Max | X_postmin - X | = "<<maxjumpr[i];
      if(maxjumpr[i]>0.02) {
        warning = true;
        std::cout<<" !! This should be zero for a properly discretized minimum energy path!\n"
        " Values greater than ~0.02 should be considered risky;"
        " Perhaps consider more integration points here. WARNING";
      }

      std::cout<<"\n\t\t        | dF_next - dF | = "<<_ddF;
      std::cout<<"\n\t\t Approx |  F_next - F  | = "<<_dF;
      std::cout<<"\n\t\t        |  r_next - r  | = "<<integr[i+1]-integr[i]<<"\n";

      if(_ddF<0.02 && std::fabs(dfer[i])<0.02) {
        warning = true;
        std::cout<<"Very flat segment! "
        "Perhaps remove knot "<<i+1<<" or "<<i+2<<"?\n. FAIL";
      }

      if(_dF>0.1) {
        warning = true;
        std::cout<<"Large free energy change! "
        "Consider more integration points here. FAIL";
      }

      ddF[0] += _ddF/double(dfer.size()-1);
      ddF[1] += _ddF*_ddF/double(dfer.size()-1);

      F += -(dfer[i]+dfer[i+1])/2. * (integr[i+1]-integr[i]);
      std::cout<<"\n -------- \n";
    }
    ddF[1] -= ddF[0]*ddF[0];
    std::cout<<"\n\t Knot "<<i+1<<": r = "<<integr[i]<<" |dF| = "<<std::fabs(dfer[i])<<" , F = "<<F<<"\n";
    std::cout<<"\n\t\t   Max | X_postmin - X | = "<<maxjumpr[i];
    std::cout<<"\n -------- \n";
    std::cout<<"\n\n\tAverage, Std in d|dF|:"<<ddF[0]<<" , "<<sqrt(ddF[1])<<std::endl;



    std::cout<<"\n\n--------------- Integration checks at zero temperature -----\n";
    std::cout<<"\n\tEnergy Barrier: "<<std::setprecision(5)<<E_bar_dense;
    std::cout<<" eV, Force Integration Barrier: "<<std::setprecision(5)<<F_bar_dense<<" eV";
    std::cout<<"\n\n\tAbsolute error of: "<<std::setprecision(5)<<std::fabs(F_bar_dense-E_bar_dense)*1000.<<" meV ";
    std::cout<<" ("<<std::setprecision(5)<<(F_bar_dense/E_bar_dense-1.0)*100.<<"%)"<<std::endl;
    if(std::fabs(F_bar_dense-E_bar_dense)<=0.005) {
      std::cout<<"\n\tError < 5meV - within likely potential accuracy, OK for sampling!\n";
    } else {
      warning=true;
      std::cout<<"\n\tError > 5meV, probably too high! FAIL\n";
    }
    if(std::fabs(F_bar_dense-E_bar_dense)>=0.001) {
      std::cout<<"\t\tTo reduce this error: ";
      std::cout<<"\n\t\t More NEB images and/or increasing nPlanes in config.xml\n\n"<<std::endl;
      std::cout<<"\n\t\t Try Rediscretize==0/1 to change how images are interpolated\n\n"<<std::endl;

    }

    if(warning) std::cout<<"\nPathway checks failed / warnings generated! See above\n"<<std::endl;
    else std::cout<<"\nPathway checks passed!\n"<<std::endl;

  }

  // close down LAMMPS instances
  sim.close();

  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
