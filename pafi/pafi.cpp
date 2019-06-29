#include "pafi.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Load input file
  Parser params("./config.xml");


  if(nProcs%params.CoresPerWorker!=0) {
    if(rank==0) std::cout<<"CoresPerWorker must factorize nProcs!\n";
    exit(-1);
  }


  // Find fresh dump folder name
  std::string new_dump_dir = params.dump_dir;
  boost::filesystem::path p(params.dump_dir);

  if(boost::filesystem::is_directory(p)) {
    for (int i = 1; boost::filesystem::is_directory(p) && i < 20; ++i) {
      std::stringstream ss;
      ss << params.dump_dir << "_" << i;
      p = ss.str();
      new_dump_dir = ss.str();
    }
    params.dump_dir = new_dump_dir;
    params.parameters["DumpFolder"] = new_dump_dir;
  }

  if(boost::filesystem::is_directory(p)) {
    if(rank==0) std::cout<<"Ran out of dump folder names!"<<std::endl;
    exit(-1);
  }

  if(rank==0) {
    boost::filesystem::create_directory(p);
    params.welcome_message();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  const int nWorkers = nProcs / params.CoresPerWorker;
  const int instance = rank / params.CoresPerWorker;
  const int local_rank = rank % params.CoresPerWorker;
  const int nRes = 7; // TODO remove this
	const int nRepeats = params.nRepeats;



  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  int *seed = new int[nWorkers];
  if(rank==0) for(int i=0;i<nWorkers;i++) seed[i]= 137*i + static_cast<int>(std::time(0))/1000;
  MPI_Bcast(seed,nWorkers,MPI_INT,0,MPI_COMM_WORLD);
  //if(local_rank==0) seed[0]=13*instance+static_cast<int>(std::time(0))/1000;
  MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Bcast(seed,1,MPI_INT,0,instance_comm);
  params.seed(seed[instance]);

  if(rank==0) std::cout<<"\n\nSet up "<<nWorkers<<" workers with   "<<params.CoresPerWorker<<" cores per worker\n\n";

  Simulator sim(instance_comm,params,rank);

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

  sim.make_path(params.KnotList);

  if(rank==0) std::cout<<"\n\nPath Loaded\n\n";

  const int vsize = 3 * sim.natoms;
  const int rsize = nRes*nWorkers*nRepeats;

  double temp, dr;
  std::string rstr,Tstr,dump_fn,fn,temp_dump_dir;
  std::ofstream raw;

  if (params.nPlanes>1) dr = (params.stopr-params.startr)/(double)(params.nPlanes-1);
  else dr = 0.1;


  for(double T = params.lowT; T <= params.highT;) {

    Tstr = boost::lexical_cast<std::string>((int)(T));

    temp_dump_dir = params.dump_dir+"/"+Tstr+"K";

    if(rank==0) {
      fn = temp_dump_dir + "/raw_ensemble_output_"+Tstr+"K";
      boost::filesystem::create_directory(temp_dump_dir);
      std::cout<<"opening dump file "<<fn<<std::endl;
      raw.open(fn.c_str(),std::ofstream::out);
      std::cout<<"\nStarting T="+Tstr+"K run\n\n";
    }
    double *local_res = new double[rsize];
    double *results = new double[nRes];
    double *local_dev = new double[vsize];
		double *local_dev_sq = new double[vsize];
    double *all_dev = NULL, *all_dev_sq = NULL, *all_res = NULL;
    std::vector<double> integr, dfer, dfere, psir;

    if (rank == 0) {
      all_dev = new double[vsize];
      all_dev_sq = new double[vsize];
      all_res = new double[rsize];
			std::cout<<"bra-kets <> == time averages,  av/err== ensemble average/error"<<std::endl;
      std::cout<<std::setw(15)<<"r";
      std::cout<<std::setw(15)<<"av(<Tpre>)";
      std::cout<<std::setw(15)<<"av(<Tpost>)";
      std::cout<<std::setw(15)<<"err(<Tpost>)";
      std::cout<<std::setw(15)<<"av(<dF/dr>)";
			std::cout<<std::setw(15)<<"err(<dF/dr>)";
      std::cout<<std::setw(15)<<"av(|<X>-U|)";
      std::cout<<std::setw(15)<<"err(|<X>-U|)";
      std::cout<<std::setw(20)<<"av(|(<X>-U).N|)";
      std::cout<<std::setw(20)<<"av(<N_true>.N)";
      std::cout<<std::setw(20)<<"err(<N_true>.N)";
      std::cout<<"\n";
    }

    for (double r = params.startr; r <= params.stopr+0.5*dr; r += dr ) {

      rstr = boost::str(boost::format("%.4f") % r);

      for(int i=0;i<rsize;i++) local_res[i] = 0.0;
			for(int i=0;i<vsize;i++) local_dev_sq[i] = 0.0;

			if(nRepeats>1 && rank==0) std::cout<<std::setw(15)<<"Repeat: "<<std::flush;
			for (int ir=0;ir<nRepeats;ir++) {
      	sim.sample(r, T, results, local_dev);
      	if(local_rank == 0) {
        	for(int i=0;i<nRes;i++) local_res[(instance*nRepeats+ir)*nRes  + i] = results[i];
					for(int i=0;i<vsize;i++) local_dev_sq[i] += local_dev[i] / double(nRepeats);
				}
				if(nRepeats>1 && rank==0) std::cout<<std::setw(15)<<ir+1<<"/"<<nRepeats<<" "<<std::flush;

			}
			if(nRepeats>1 && rank==0) std::cout<<std::endl;
			for(int i=0;i<vsize;i++) local_dev[i] = local_dev_sq[i];
			for(int i=0;i<vsize;i++) local_dev_sq[i] = local_dev[i]*local_dev[i];


      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(local_res,all_res,rsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      MPI_Reduce(local_dev,all_dev,vsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Reduce(local_dev_sq,all_dev_sq,vsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

      if (rank==0) {

        double *final_res = new double[2*nRes];
        dump_fn = temp_dump_dir+"/dev_"+rstr+"_"+Tstr+"K.dat";

        sim.write_dev(dump_fn,r,all_dev,all_dev_sq);

				// (instance*nRepeats+ir)*nRes  + i

        for(int i=0;i<nWorkers*nRepeats;i++) raw<<all_res[i*nRes+2]<<" ";
        raw<<std::endl;

        for(int j=0;j<2*nRes;j++) final_res[j] = 0.;

        for(int i=0;i<nWorkers*nRepeats;i++) for(int j=0;j<nRes;j++)
          final_res[j] += all_res[i*nRes+j] / double(nWorkers*nRepeats);

        for(int i=0;i<nWorkers*nRepeats;i++) for(int j=0;j<nRes;j++) {
          temp = all_res[i*nRes+j]-final_res[j];
          final_res[j+nRes] += temp * temp / double(nWorkers*nRepeats);
        }
				// under assumption that sample time is longer than sample autocorrelations,
				// expected errors in time averages is ensemble average variance divided by nWorkers
				// raw_output gives data to confirm this assumption (CLT with grouping)
        for(int j=0;j<nRes;j++) final_res[j+nRes] = sqrt(final_res[j+nRes]/double(nWorkers*nRepeats));

        integr.push_back(r);
        dfer.push_back(final_res[2]); // <dF/dr>


				dfere.push_back(final_res[2+nRes]);//"err(<dF/dr>)"
        psir.push_back(final_res[4]); // <Psi>
        std::cout<<std::setw(15)<<r;//"r"
        std::cout<<std::setw(15)<<final_res[0];//"av(<Tpre>)"
        std::cout<<std::setw(15)<<final_res[1];//"av(<Tpost>)"
        std::cout<<std::setw(15)<<final_res[1+nRes];//"std(<Tpost>)"
        std::cout<<std::setw(15)<<final_res[2];//"av(<dF/dr>)"
				std::cout<<std::setw(15)<<final_res[2+nRes];//"err(<dF/dr>)"
        std::cout<<std::setw(15)<<final_res[6];//"av(|<X>-U|)"
        std::cout<<std::setw(15)<<final_res[6+nRes];//"err(|<X>-U|)"
        std::cout<<std::setw(20)<<final_res[5];//"av(|<X>-U).(dU/dr)|)"
        std::cout<<std::setw(20)<<final_res[4];//"av(Psi)"
        std::cout<<std::setw(20)<<final_res[4+nRes];//"std(Psi)"
        std::cout<<"\n";
      }
    }

    // free_energy_profile
    if(rank==0){
      raw.close();
      std::cout<<"T="+Tstr+"K run complete, integrating FEP....\n\n";
      spline dfspl,psispl,dfespl;

      dfspl.set_points(integr,dfer);
      dfespl.set_points(integr,dfere);
      psispl.set_points(integr,psir);

      std::vector<std::array<double,4>> fF;
      std::array<double,4> fline;
      double dr = 1./(double)(30 * integr.size());
      double rtF=0.0;

      fline[0]=params.startr; fline[1]=rtF; fline[2]=0.; fline[3] = psispl(params.startr);
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
      fn = temp_dump_dir + "/free_energy_profile_"+Tstr+"K";
      out.open(fn.c_str(),std::ofstream::out);
      out<<"# r F(r) av(<dF/dr>) err(<dF/dr>) av(<Psi>)\n";
      for(auto l: fF)
        out<<l[0]<<" "<<l[1]+l[2]<<" "<<dfspl(l[0])<<" "<<dfespl(l[0])<<" "<<l[3]<<"\n";
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
