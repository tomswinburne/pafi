#include "pafi.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Load input file
  Parser params("./config.xml",false);

  if(!params.xml_success) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Configuration file could not be read!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }

  if(nProcs%params.CoresPerWorker!=0) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"CoresPerWorker must factorize nProcs!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }



  // Find fresh dump folder name - no nice solution here as
  // directory creation requires platform dependent features
  // which we omit for portability
  int *int_dump_suffix = new int[1];
  std::ofstream raw;
  std::string params_file =
    params.dump_dir+"/params_"+std::to_string(int_dump_suffix[0]);

  // try to write to a file with a unique suffix
  if(rank==0) {
    for (int_dump_suffix[0]=0; int_dump_suffix[0] < 100; int_dump_suffix[0]++) {
      params_file =
        params.dump_dir+"/params_"+std::to_string(int_dump_suffix[0]);
      std::cout<<"Testing for existence of "<<params_file<<".... ";
      if(!file_exists(params_file)) {
        std::cout<<"it doesn't!, trying to open for writing.... ";
        raw.open(params_file.c_str(),std::ofstream::out);
        if(raw.is_open()) {
          std::cout<<" done!\n";
          raw<<params.welcome_message();
          raw.close();
          break;
        } else std::cout<<"cannot!\n";
      } else std::cout<<"it exists!\n";
    }
  }
  MPI_Bcast(int_dump_suffix,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if(int_dump_suffix[0]==100) {
    if(rank==0) {
      std::cout<<"\n\n\n*****************************\n\n\n";
      std::cout<<"Could not write to output path / find directory! Exiting!\n";
      std::cout<<"\n\n\n*****************************\n\n\n";
    }
    exit(-1);
  }

  if(rank==0)  std::cout<<params.welcome_message();

  MPI_Barrier(MPI_COMM_WORLD);

  const int nWorkers = nProcs / params.CoresPerWorker;
  const int instance = rank / params.CoresPerWorker;
  const int local_rank = rank % params.CoresPerWorker;
  const int nRes = 8; // TODO remove this
  const int nRepeats = params.nRepeats;



  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  int *seed = new int[nWorkers];
  if(rank==0) for(int i=0;i<nWorkers;i++)
    seed[i]= 137*i + static_cast<int>(std::time(0))/1000;
  MPI_Bcast(seed,nWorkers,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  params.seed(seed[instance]);

  if(rank==0) std::cout<<"\n\nInitializing "<<nWorkers<<" workers "
                      "with "<<params.CoresPerWorker<<" cores per worker\n\n";

  Simulator sim(instance_comm,params,instance,nRes);

  if(!sim.has_pafi) {
    if(rank==0)
      std::cout<<"PAFI Error: missing "<<sim.pafi_package<<" package in LAMMPS"<<std::endl;
    exit(-1);
  }
  if(sim.error_count>0 && local_rank==0) std::cout<<sim.last_error()<<std::endl;

  if (rank == 0) {
    std::cout<<"Loaded input data of "<<sim.natoms<<" atoms\n";
    std::cout<<"Supercell Matrix:\n";
    auto cell = sim.getCellData();
    for(int i=0;i<3;i++) {
      std::cout<<"\t";
      for(int j=0;j<3;j++) std::cout<<sim.pbc.cell[i][j]<<" ";
      std::cout<<"\n";
    }
    std::cout<<"\n\n";
  }

  sim.make_path(params.PathwayConfigurations,params.real_coord);
  if(sim.error_count>0 && local_rank==0) std::cout<<sim.last_error()<<std::endl;

  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";

  const int vsize = 3 * sim.natoms;
  const int rsize = nRes*nWorkers;
  bool not_finished_sampling;
  int total_valid_data, totalRepeats, total_invalid_data;
  double temp, t_max_jump, p_jump, f_mean;
  double *f_error = new double[1];
  std::string rstr,Tstr,dump_fn,fn,dump_suffix;
  std::map<std::string,double> results;

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
  std::vector<double> integr, dfer, dfere, psir;
  std::vector<double> valid_res, invalid_res;
  std::list<int> raw_dump_indicies = {2,3,4,7};// <f> std(f) <Psi> <Jump>


  std::vector<double> sample_r = params.sample_r(sim.pathway_r);
  

  for(double T = params.lowT; T <= params.highT;) {

    Tstr = std::to_string((int)(T));

    dump_suffix = "_"+Tstr+"K_"+std::to_string(int_dump_suffix[0]);

    if(rank==0) {
      fn = params.dump_dir + "/raw_ensemble_output"+dump_suffix;
      std::cout<<"opening dump file "<<fn<<std::endl;
      raw.open(fn.c_str(),std::ofstream::out);
      std::cout<<"\nStarting T="+Tstr+"K run\n\n";
      if(std::stoi(params.parameters["OverDamped"])==1) {
        std::cout<<"CAUTION : when OverDamped==1 Tpre/post are estimated from "
        "equipartition. This will be less accurate at high temperature.\n\n";
      }
    }

    for(int i=0;i<nWorkers;i++) valid[i]=0;
    for(int i=0;i<rsize;i++) local_res[i] = 0.0;
    for(int i=0;i<vsize;i++) local_dev[i] = 0.0;
    for(int i=0;i<vsize;i++) local_dev_sq[i] = 0.0;
    results.clear();
    integr.clear();
    dfer.clear();
    dfere.clear();
    psir.clear();
    valid_res.clear();
    invalid_res.clear();

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


    for(auto r: sample_r) {
    //for (double r = params.startr; r <= params.stopr+0.5*dr; r += dr ) {
      valid_res.clear();
      invalid_res.clear();
      rstr = std::to_string(r);
      for(int i=0;i<rsize;i++) local_res[i] = 0.0;
      for(int i=0;i<vsize;i++) local_dev[i] = 0.0;
      for(int i=0;i<vsize;i++) local_dev_sq[i] = 0.0;


      total_valid_data=0;
      totalRepeats=0;
      t_max_jump=0.0;
      not_finished_sampling = true;
      while (not_finished_sampling) {

        sim.sample(r, T, results, local_dev);
        double aa=0.0;

        for(int i=0;i<vsize;i++) local_dev_sq[i] = local_dev[i]*local_dev[i];
        for(int i=0;i<vsize;i++) aa=std::max(aa,local_dev_sq[i]);
        totalRepeats++;

        if(local_rank == 0) {
          if(sim.error_count>0) std::cout<<sim.last_error()<<std::endl;
          local_res[instance*nRes + 0] = results["preT"];
          local_res[instance*nRes + 1] = results["postT"];
          local_res[instance*nRes + 2] = results["aveF"];
          local_res[instance*nRes + 3] = results["stdF"];
          local_res[instance*nRes + 4] = results["avePsi"];
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
        if(valid[instance]==0 && !params.postMin)
          for(int i=0;i<vsize;i++) {
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
          if(rank==0) std::cout<<"Reference path too unstable for sampling.\n"
                                "Try a lower temperature. See README for tips";
          break;
        }
        not_finished_sampling = bool(total_valid_data<=int(params.redo_thresh*nWorkers*nRepeats));

        // determine force error
        f_error[0]=0.0;
        if(rank==0){
          f_mean = 0.0;
          for (int i=0;i<total_valid_data;i++)
            f_mean += valid_res[i*nRes+2] / double(total_valid_data);
          for (int i=0;i<total_valid_data;i++) {
            // we have N^2 as require standard error on mean
            temp = (f_mean-valid_res[i*nRes+2]) / double(total_valid_data);
            f_error[0] += temp * temp ;
          }
          f_error[0] = std::sqrt(f_error[0]);
        }
        MPI_Bcast(f_error,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        not_finished_sampling += bool(f_error[0] >= params.f_error_thresh);

        if(rank==0) std::cout<<"#"<<std::flush;
      }
      // collate
      if(rank==0) {
        total_invalid_data = totalRepeats*nWorkers - total_valid_data;

        // deviation vectors
        if(params.postMin) for(int i=0;i<vsize;i++) {
          all_dev[i]/=double(totalRepeats*nProcs);
          all_dev_sq[i]/=double(totalRepeats*nProcs);
        } else if(total_valid_data>0) for(int i=0;i<vsize;i++) {
          all_dev[i]/=double(total_valid_data*params.CoresPerWorker);
          all_dev_sq[i]/=double(total_valid_data*params.CoresPerWorker);
        }

        dump_fn = params.dump_dir+"/dev_"+rstr+dump_suffix+".dat";
        sim.write_dev(dump_fn,r,all_dev,all_dev_sq);

        std::cout<<" "<<totalRepeats<<"/"<<nRepeats<<" ";
        std::cout<<total_valid_data<<"/"<<nRepeats*nWorkers<<" ";
        std::cout<<valid_res.size()/nWorkers<<std::endl;
        // raw output
        raw<<r<<" "<<total_valid_data<<" ";
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
          final_res[j+nRes] += temp * temp / double(total_valid_data);
        }
        // under assumption that sample time is longer than autocorrelations,
        // expected variance in time average is ensemble variance / nWorkers
        // raw_output gives data to confirm this assumption (CLT with grouping)
        if(total_valid_data>0) for(int j=0;j<nRes;j++)
          final_res[j+nRes] = sqrt(final_res[j+nRes])/sqrt(double(total_valid_data));

        integr.push_back(r);
        dfer.push_back(final_res[2]); // <dF/dr>
        dfere.push_back(final_res[2+nRes]);//"err(<dF/dr>)"
        psir.push_back(final_res[4]); // <Psi>
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
    if(rank==0){
      raw.close();
      std::cout<<"\nT="+Tstr+"K run complete, integrating FEP....\n\n";
      spline dfspl,psispl,dfespl;

      dfspl.set_points(integr,dfer);
      dfespl.set_points(integr,dfere);
      psispl.set_points(integr,psir);

      std::vector<std::array<double,4>> fF;
      std::array<double,4> fline;
      double dr = 1./(double)(2 * integr.size());
      double rtF=0.0;

      fline[0]=params.startr;
      fline[1]=rtF;
      fline[2]=0.;
      fline[3] = psispl(params.startr);
      fF.push_back(fline);

      for(double sr=params.startr+dr; sr<= params.stopr+dr*0.5; sr+=dr) {
        rtF -=  dr * dfspl(sr);
        fline[0] = sr;
        fline[1] = rtF;
        fline[2] = BOLTZ*T*log(fabs(psispl(sr)/psispl(0.)));
        fline[3] = psispl(sr);
        fF.push_back(fline);
      }
      std::ofstream out;
      fn = params.dump_dir + "/free_energy_profile"+dump_suffix;
      out.open(fn.c_str(),std::ofstream::out);
      out<<"# r F(r) av(<dF/dr>) err(<dF/dr>) av(<Psi>)\n";
      for(auto l: fF) {
        out<<l[0]<<" "<<l[1]<<" ";
        out<<dfspl(l[0])<<" "<<dfespl(l[0])<<" "<<l[3]<<"\n";
      }
      out.close();
      std::cout<<"Integration complete\n\n";
    }

    if( params.highT > params.lowT + 1.0 && params.TSteps > 1) \
      T += ( params.highT - params.lowT ) / ( params.TSteps - 1 );
    else break;
  }

  // close down LAMMPS instances
  sim.close();

  // close down MPI
  MPI_Comm_free(&instance_comm);

  MPI_Finalize();

}
