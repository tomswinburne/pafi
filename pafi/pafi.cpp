#include "pafi.hpp"

int main(int narg, char **arg) {
  MPI_Init(&narg,&arg);
  int rank, nProcs;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nProcs);

  // Parse input file
  Parser params("./config.xml");

  double T = boost::lexical_cast<double>(params.parameters["LowTemperature"]);


  if(rank==0) std::cout << "PAFI: XML read\n";

  if(nProcs%params.CoresPerWorker!=0) {
    if(rank==0) std::cout<<"CoresPerWorker must factorize nProcs!\n";
    exit(-1);
  }

  const int nWorkers = nProcs / params.CoresPerWorker;
  const int instance = rank / params.CoresPerWorker;
  const int local_rank = rank % params.CoresPerWorker;
  const int nRes = 7;



  MPI_Comm instance_comm;
  MPI_Comm_split(MPI_COMM_WORLD,instance,0,&instance_comm);

  if(rank==0) std::cout<<"\nSet up "<<nWorkers<<" workers with   "<<params.CoresPerWorker<<" cores per worker\n";

  Simulator sim(instance_comm,params,rank);

  if (rank == 0) {
    std::cout<<"Loaded input data of "<<sim.natoms<<"atoms\n";
    std::cout<<"Supercell Matrix:\n";
    auto cell = sim.getCellData();
    for(int i=0;i<3;i++) {
      std::cout<<"\t";
      for(int j=0;j<3;j++) std::cout<<sim.pbc.cell[i][j]<<" ";
      std::cout<<"\n";
    }
  }

  sim.make_path(params.KnotList);

  if(rank==0) std::cout<<"\n\nPATH LOADED\n";

  const int vsize = 3 * sim.natoms;
  const int rsize = nRes*nWorkers;


  double *local_res = new double[rsize];
  double *results = new double[nRes];
  double *local_dev = new double[vsize];
  double *all_dev = NULL, *all_dev_sq = NULL, *all_res = NULL;
  std::vector<double> integr, dfer, dfere, psir;

  if (rank == 0) {
    all_dev = new double[vsize];
    all_dev_sq = new double[vsize];
    all_res = new double[rsize];

  }

  for (double r=0.00; r<=1.; r+=0.1 ) {

    for (int i=0;i<rsize;i++) local_res[i] = 0.;

    sim.sample(r, T, results, local_dev);

    if(local_rank == 0)
      for(int i=0;i<nRes;i++) local_res[instance*nRes + i] = results[i];

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(local_res,all_res,rsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Reduce(local_dev,all_dev,vsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<vsize;i++) local_dev[i] *= local_dev[i];
    MPI_Reduce(local_dev,all_dev_sq,vsize,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if (rank==0) {
      std::string dfn = "dumps/dev_"+boost::lexical_cast<std::string>(r)+".dat";
      sim.write_dev(dfn,r,all_dev,all_dev_sq);

      double *final_res = new double[2*nRes];

      for(int j=0;j<2*nRes;j++) final_res[j] = 0.;

      for(int i=0;i<nWorkers;i++) for(int j=0;j<nRes;j++)
        final_res[j] += all_res[i*nRes+j] / (double)nWorkers;
      double temp;
      for(int i=0;i<nWorkers;i++) for(int j=0;j<nRes;j++) {
        temp = all_res[i*nRes+j]-final_res[j];
        final_res[j+nRes] += temp * temp / (double)nWorkers;
      }
      for(int j=0;j<nRes;j++) final_res[j+nRes] = sqrt(final_res[j+nRes]);

      integr.push_back(r);
      dfer.push_back(final_res[2]); // <dF/dr>
      dfere.push_back(final_res[2+nRes]); // std(dF/dr)
      psir.push_back(final_res[4]); // <Psi>

      std::cout<<"SAMPLED: r:"<<r<<" T:"<<T<<"\n";

      std::cout<<"\t <T>: "<<final_res[1]<<" +/- "<<final_res[1+nRes]<<"\n";
      std::cout<<"\t <dF/dr>: "<<final_res[2]<<" +/- "<<final_res[2+nRes]<<"\n";
      std::cout<<"\t <|dX|>: "<<final_res[6]<<" +/- "<<final_res[6+nRes]<<"\n";
      std::cout<<"\n\n\n";
    }
  }

  // close down LAMMPS instances
  sim.close();


  // Initial dump file,
  if(rank==0){
    std::cout<<"Simulation complete, integrating free energy profile....\n\n";
    spline dfspl,psispl,dfespl;

    dfspl.set_points(integr,dfer);
    dfespl.set_points(integr,dfere);
    psispl.set_points(integr,psir);

    std::vector<std::array<double,3>> fF;
    std::array<double,3> fline;
    double dr = 1./(double)(30 * integr.size());
    double rtF=0.0;

    fline[0]=0.; fline[1]=rtF; fline[2]=0.;
    fF.push_back(fline);
    for(double sr=dr; sr<=1.0; sr+=dr) {
      rtF -=  dr * dfspl(sr);
      fline[0]=sr; fline[1]=rtF; fline[2]=BOLTZ*T*log(psispl(sr)/psispl(0.));
      fF.push_back(fline);
    }
    std::ofstream out;
    std::string fn = "free_energy_profile";
    out.open(fn.c_str(),std::ofstream::out);
    out<<"# r F(r) <dF/dr> std(dF/dr)\n";
    for(auto l: fF)
      out<<l[0]<<" "<<l[1]<<" "<<l[2]<<" "<<dfspl(l[0])<<" "<<dfespl(l[0])<<"\n";
    out.close();
  }

  // close down MPI
  MPI_Comm_free(&instance_comm);



  MPI_Finalize();

}
