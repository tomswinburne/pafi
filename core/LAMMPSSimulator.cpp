#include "LAMMPSSimulator.hpp"

LAMMPSSimulator::LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int rank) : GeneralSimulator::GeneralSimulator(instance_comm,p,rank) {
  int i;
  tag = rank;
  params = &p;
  comm = &instance_comm;

  // set up LAMMPS
  std::string logfile = params->loglammps ? ("log.lammps."+std::to_string(tag)) : "none";
  char **lmparg = new char*[5];
  lmparg[0] = NULL; // required placeholder for program name
  lmparg[1] = (char *) "-screen";
  lmparg[2] = (char *) "none";
  lmparg[3] = (char *) "-log";
  lmparg[4] = (char *) logfile.c_str();
  lmp = new LAMMPS_NS::LAMMPS(5,lmparg,*comm);

  run_script("Input");
  natoms = *((int *) lammps_extract_global(lmp,(char *) "natoms"));

  nlocal = static_cast<int> (lmp->atom->nlocal);
  natoms = static_cast<int> (lmp->atom->natoms);

  // for local splining info. Other ways obviously, but this could have less comm?



  int *lid = (int *)lammps_extract_atom(lmp,(char *) "id");
  for(i=0;i<nlocal;i++) loc_id.push_back(lid[i]-1);
  for(i=0;i<3*nlocal;i++) loc_data.push_back(0.0);
  for(i=0;i<3*natoms;i++) glo_data.push_back(0.0);

  // get cell info
  pbc.load(getCellData());

  /*
  species = std::vector<int>(natoms,1);
  lammps_gather_atoms(lmp,(char *) "type",0,1,&species[0]);
  image = std::vector<int>(natoms,0);
  lammps_gather_atoms(lmp,(char *) "image",0,1,&image[0]);
  id = std::vector<int>(natoms,0);
  lammps_gather_atoms(lmp,(char *) "id",0,1,&id[0]);
  q = std::vector<double>(natoms,0.0);
  lammps_gather_atoms(lmp,(char *) "q",0,1,&q[0]);
  */
  delete [] lmparg;

};

/*
  Load xyz configuration from data file and put in vector
*/
void LAMMPSSimulator::load_config(std::string file_string,std::vector<double> &x) {
  lammps_command(lmp,(char *)"delete_atoms group all");
  std::string fn = "read_data ";
  fn += file_string;
  fn += " add merge";
  lammps_command(lmp,(char *)fn.c_str());
  lammps_gather_atoms(lmp,(char *) "x",1,3,&x[0]);
};

/*
  Parse and run script from configuration file
*/
void LAMMPSSimulator::run_script(std::string sn){
  std::vector<std::string> strv = params->Script(sn);
  for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
};

/*
  Parse and run script from string with linebreaks
*/
void LAMMPSSimulator::run_commands(std::vector<std::string> strv) {
  for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
};

void LAMMPSSimulator::run_commands(std::string strv) {
  run_commands(params->Parse(strv));
};


/*
  Fill configuration, path, tangent and tangent gradient. Return tangent norm
*/
void LAMMPSSimulator::populate(double r, double *scale, double &nm, bool justpos) {
  int i,j;
  position=r;
  std::vector<double> t(3*natoms,0.);
  for(j=0;j<3;j++) norm_com[j] = 0.0;
  nm=0.;

  GeneralSimulator::fill_path(r,0,t);
  for(i=0;i<3*natoms;i++) t[i] *= scale[i%3];
  lammps_scatter_atoms(lmp,(char *)"x",1,3,&t[0]);

  if(!justpos) {
    lammps_scatter_atoms(lmp,(char *)"path",1,3,&t[0]);

    fill_path(r,1,t,scale);

    // Center of mass projection and tangent length
    for(i=0;i<3*natoms;i++) norm_com[i%3] += t[i]/(double)natoms;
    for(i=0;i<3*natoms;i++) t[i] -= norm_com[i%3];
    for(i=0;i<3*natoms;i++) nm += t[i]*t[i];
    nm = sqrt(nm);

    for(i=0;i<3*natoms;i++) t[i] /= nm;
    lammps_scatter_atoms(lmp,(char *)"norm",1,3,&t[0]);


    fill_path(r,2,t,scale);
    for(i=0;i<3*natoms;i++) t[i] /= nm * nm;
    lammps_scatter_atoms(lmp,(char *)"dnorm",1,3,&t[0]);

  }

  lammps_command(lmp,(char *)"run 0");
};

void LAMMPSSimulator::fill_path(double r,int der,std::vector<double> &vec, double *scale) {
  GeneralSimulator::fill_path(r,der,vec);
  for(int i=0;i<vec.size();i++) vec[i] *= scale[i%3];
};

/*
  Rescale simulation cell
*/
void LAMMPSSimulator::rescale_cell(double *scale) {
  std::string cmd, ssx,ssy,ssz;
  ssx = std::to_string(scale[0]);
  ssy = std::to_string(scale[1]);
  ssz = std::to_string(scale[2]);
  cmd ="change_box all x scale "+ssx+" y scale "+ssy+" z scale "+ssz+"\n";
  cmd += "run 0";
  run_commands(cmd);
};

/*
  Main sample run. Results vector should have thermalization temperature,
  sample temperature <f>, <f^2>, <psi> and <x-u>.n
*/
void LAMMPSSimulator::sample(double r, double T, double *results, double *dev) {
  position = r;
  temperature = T;
  int i,j;
  std::string cmd;
  double sampleT, dm;
  double *lmp_ptr;
  double **lmp_dev_ptr;

  double scale[3] = {1.0,1.0,1.0};

  populate(r,scale,norm_mag,true);


  // Stress Fixes
  run_script("PreRun");
  cmd = "run 0";
  run_commands(cmd);

  expansion(T,scale);
  rescale_cell(scale);

  populate(r,scale,norm_mag);


  params->parameters["Temperature"] = std::to_string(T);
  cmd = "fix hp all hp "+params->parameters["Temperature"]+" ";
  cmd += params->parameters["Friction"]+" ";
  cmd += params->seed_str()+" overdamped ";
  cmd += params->parameters["OverDamped"]+" com 1\nrun 0";
  run_commands(cmd);

  refE = getEnergy();
  cmd = "reset_timestep 0\n";
  cmd += "fix ae all ave/time 1 "+params->parameters["ThermWindow"]+" ";
  cmd += params->parameters["ThermSteps"]+" v_pe\n";
  cmd += "run "+params->parameters["ThermSteps"];
  run_commands(cmd);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
  sampleT = (*lmp_ptr-refE)/natoms/1.5/8.617e-5;
  cmd = "unfix ae\nrun 0";
  run_commands(cmd);
  results[0] = sampleT;

  std::string SampleSteps = params->parameters["SampleSteps"];
  cmd = "reset_timestep 0\n";
  cmd += "fix ae all ave/time 1 "+SampleSteps+" "+SampleSteps+" v_pe\n";
  cmd += "fix ad all ave/deviation 1 "+SampleSteps+" "+SampleSteps+"\n";
  cmd += "fix af all ave/time 1 "+SampleSteps+" "+SampleSteps+" f_hp[1] f_hp[2] f_hp[3] f_hp[4]\n";
  cmd += "run "+SampleSteps;
  run_commands(cmd);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
  sampleT = (*lmp_ptr-refE)/natoms/1.5/8.617e-5;
  results[1] = sampleT;
  lammps_free(lmp_ptr);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,0,0);
  results[2] = *lmp_ptr * norm_mag;
  lammps_free(lmp_ptr);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,1,0);
  results[3] = *lmp_ptr * norm_mag * norm_mag - results[2] * results[2];
  lammps_free(lmp_ptr);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,2,0);
  results[4] = *lmp_ptr;
  lammps_free(lmp_ptr);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,3,0);
  results[5] = *lmp_ptr;
  lammps_free(lmp_ptr);

  lammps_gather_peratom_fix(lmp,(char *)"ad",3,dev);
  dm=0.; for(i=0;i<3*natoms;i++) dm += dev[i] * dev[i];
  results[6] = dm;

  cmd = "unfix ae\nunfix af\nunfix ad\nunfix hp";
  run_commands(cmd);

  // Stress Fixes
  run_script("PostRun");
  cmd = "run 0"; run_commands(cmd);

  // rescale back
  for(j=0;j<3;j++) scale[j] = 1.0/scale[j];

  populate(r,scale,norm_mag);
  rescale_cell(scale);

};

double LAMMPSSimulator::getEnergy() {
  lammps_command(lmp,(char *) "run 0");
  double * lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  if (lmpE == NULL) {
    lammps_command(lmp,(char *) "compute pe all pe");
    lammps_command(lmp,(char *) "variable pe equal pe");
    lammps_command(lmp,(char *) "run 0");
    lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  }
  double baseE = *lmpE;
  return baseE;
};

double LAMMPSSimulator::getForceEnergy(double *f) {
  lammps_command(lmp,(char *) "run 0");
  double * lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  if (lmpE == NULL) {
    lammps_command(lmp,(char *) "compute pe all pe");
    lammps_command(lmp,(char *) "variable pe equal pe");
    lammps_command(lmp,(char *) "run 0");
    lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
  }
  lammps_gather_atoms(lmp,(char *) "f",1,3,f);
  double baseE = *lmpE;
  return baseE;
};

// Fill 9D array with Lx, Ly, Lz, xy, xz, yz, then periodicity in x, y, z
std::array<double,9> LAMMPSSimulator::getCellData() {
  int i;
  std::array<double,9> cell;
  for(i=0; i<9; i++) cell[i] = 0.;
  double *lmp_ptr;
  double *boxlo = (double *) lammps_extract_global(lmp,(char *) "boxlo");
  double *boxhi = (double *) lammps_extract_global(lmp,(char *) "boxhi");
  int *pv = (int *) lammps_extract_global(lmp,(char *) "periodicity");

  std::array<std::string,3> od;
  od[0]="xy"; od[1]="xz"; od[2]="yz";
  for(i=0; i<3; i++) {
    cell[i] = boxhi[i]-boxlo[i];
    lmp_ptr = (double *) lammps_extract_global(lmp,(char *)od[i].c_str());
    cell[3+i] = *lmp_ptr;
    cell[6+i] = pv[i];
  }
  return cell;
};

void LAMMPSSimulator::close() {
  // close down LAMMPS
  delete lmp;
};


void LAMMPSSimulator::write_dev(std::string fn, double *dev, double *dev_sq){
  int i,j;
  double scale[3];
  expansion(temperature,scale);
  std::cout<<"write_dev "<<local_rank<<" "<<position<<" "<<temperature<<" "<<scale[0]<<" "<<scale[1]<<" "<<scale[2]<<std::endl;
  std::vector<double> t(3*natoms,0.);
  fill_path(position,0,t,scale);

  if(local_rank==0) {
    std::ofstream out;
    out.open(fn.c_str(),std::ofstream::out);

    out<<"# PAFI DUMP FILE. Reference path u(r) is a Nx3 vector.\n";
    out<<"# For i=0,1,2: u_i(r) , < x_i-u_i | r > , <(x_i-u_i)^2 | r >)\n";
    for(i=0;i<natoms;i++) {
      out<<i<<" ";
      for(j=0;j<3;j++) out<<t[3*i+j]<<" ";
      for(j=0;j<3;j++) out<<dev[3*i+j]<<" ";
      for(j=0;j<3;j++) out<<dev_sq[3*i+j]<<" ";
      out<<std::endl;
    }
    out.close();
  }
};
