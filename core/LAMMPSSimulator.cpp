#include "LAMMPSSimulator.hpp"

LAMMPSSimulator::LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int rank) {
  tag = rank;
  params = &p;

  // set up LAMMPS
  char str1[32];
  char **lmparg = new char*[5];
  lmparg[0] = NULL; // required placeholder for program name
  lmparg[1] = (char *) "-screen";
  lmparg[2] = (char *) "none";
  lmparg[3] = (char *) "-log";
  //lmparg[4] = (char *) "none";
  sprintf(str1,"log.lammps.%d",rank);
  lmparg[4] = str1;

  lammps_open(5,lmparg,instance_comm,(void **) &lmp);
  run_script("Input");
  natoms = *((int *) lammps_extract_global(lmp,(char *) "natoms"));

  id = std::vector<int>(natoms,0);
  lammps_gather_atoms(lmp,(char *) "id",0,1,&id[0]);

  // get cell info
  pbc.load(getCellData());

  // Get type / image info
  species = std::vector<int>(natoms,1);
  lammps_gather_atoms(lmp,(char *) "type",0,1,&species[0]);
  s_flag=true;

  image = std::vector<int>(natoms,1);
  lammps_gather_atoms(lmp,(char *) "image",0,1,&image[0]);

  // prepare for fix?

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

/*
  Fill configuration, path, tangent and tangent gradient. Return tangent norm
*/
void LAMMPSSimulator::populate(double r, double scale, double &norm_mag) {
  std::vector<double> t(3*natoms,0.);
  double ncom[]={0.,0.,0.};
  norm_mag=0.;

  for(int i=0;i<3*natoms;i++) t[i] = pathway[i](r) * scale;
  lammps_scatter_atoms(lmp,(char *)"x",1,3,&t[0]);
  lammps_scatter_atoms(lmp,(char *)"path",1,3,&t[0]);

  for(int i=0;i<3*natoms;i++) t[i] = pathway[i].deriv(1,r) * scale;

  // Center of mass projection and tangent length
  for(int i=0;i<3*natoms;i++) ncom[i%3] += t[i]/(double)natoms;
  for(int i=0;i<3*natoms;i++) t[i] -= ncom[i%3];
  for(int i=0;i<3*natoms;i++) norm_mag += t[i]*t[i];
  norm_mag = sqrt(norm_mag);

  for(int i=0;i<3*natoms;i++) t[i] /= norm_mag;
  lammps_scatter_atoms(lmp,(char *)"norm",1,3,&t[0]);

  for(int i=0;i<3*natoms;i++) \
    t[i] = pathway[i].deriv(2,r) * scale / norm_mag / norm_mag;
  lammps_scatter_atoms(lmp,(char *)"dnorm",1,3,&t[0]);
  lammps_command(lmp,(char *)"run 0");
};

/*
  Rescale simulation cell
*/
void LAMMPSSimulator::rescale_cell(double scale) {
  std::string cmd, ss;
  ss = boost::lexical_cast<std::string>(scale);
  cmd ="change_box all x scale "+ss+" y scale "+ss+" z scale "+ss+"\n";
  cmd += "run 0";
  run_commands(params->Parse(cmd));
};

/*
  Main sample run. Results vector should have thermalization temperature,
  sample temperature <f>, <f^2>, <psi> and <x-u>.n
*/
void LAMMPSSimulator::sample(double r, double T, double *results, double *dev) {
  std::string cmd;
  double norm_mag, scale, sampleT, dm;
  double *lmp_ptr;
  double **lmp_dev_ptr;

  scale = expansion(T);
  rescale_cell(scale);
  populate(r,scale,norm_mag);

  // Stress Fixes
  run_script("PreRun");
  cmd = "run 0"; run_commands(params->Parse(cmd));


  params->parameters["Temperature"] = boost::lexical_cast<std::string>(T);
  cmd = "fix hp all hp %Temperature% %Friction% %RANDOM% overdamped ";
  cmd += params->parameters["OverDamped"]+" com 1\nrun 0";
  run_commands(params->Parse(cmd));

  refE = getEnergy();
  cmd = "reset_timestep 0\n";
  cmd += "fix ae all ave/time 1 %ThermWindow% %ThermSteps% v_pe\n";
  cmd += "run %ThermSteps%";
  run_commands(params->Parse(cmd));

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
  sampleT = (*lmp_ptr-refE)/natoms/1.5/8.617e-5;
  cmd = "unfix ae\nrun 0";
  run_commands(params->Parse(cmd));
  results[0] = sampleT;

  cmd = "reset_timestep 0\n";
  cmd += "fix ae all ave/time 1 %SampleSteps% %SampleSteps% v_pe\n";
  cmd += "fix ad all ave/deviation 1 %SampleSteps% %SampleSteps%\n";
  cmd += "fix af all ave/time 1 %SampleSteps% %SampleSteps% f_hp[1] f_hp[2] f_hp[3] f_hp[4]\nrun %SampleSteps%";
  run_commands(params->Parse(cmd));

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
  dm=0.; for(int i=0;i<3*natoms;i++) dm += dev[i] * dev[i];
  results[6] = dm;

  cmd = "unfix ae\nunfix af\nunfix ad\nunfix hp";
  run_commands(params->Parse(cmd));

  // Stress Fixes
  run_script("PostRun");
  cmd = "run 0"; run_commands(params->Parse(cmd));

  // rescale back
  scale = 1. / scale;
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

// Fill 9D array with Lx, Ly, Lz, xy, xz, yz, then periodicity in x, y, z
std::array<double,9> LAMMPSSimulator::getCellData() {
  std::array<double,9> cell;
  for(int i=0; i<9; i++) cell[i] = 0.;
  double *lmp_ptr;
  double *boxlo = (double *) lammps_extract_global(lmp,(char *) "boxlo");
  double *boxhi = (double *) lammps_extract_global(lmp,(char *) "boxhi");
  int *pv = (int *) lammps_extract_global(lmp,(char *) "periodicity");

  std::array<std::string,3> od;
  od[0]="xy"; od[1]="xz"; od[2]="yz";
  for(int i=0; i<3; i++) {
    cell[i] = boxhi[i]-boxlo[i];
    lmp_ptr = (double *) lammps_extract_global(lmp,(char *)od[i].c_str());
    cell[3+i] = *lmp_ptr;
    cell[6+i] = pv[i];
  }
  return cell;
};

void LAMMPSSimulator::close() {
  lammps_close(lmp);
};


std::string LAMMPSSimulator::header() {
  std::string res="LAMMPS dump file\n\n";
  res += boost::lexical_cast<std::string>(natoms)+" atoms\n";
  res += "1 atom types\n\n"; // TODO species
  res += "0. "+boost::lexical_cast<std::string>(pbc.cell[0][0])+" xlo xhi\n";
  res += "0. "+boost::lexical_cast<std::string>(pbc.cell[1][1])+" ylo yhi\n";
  res += "0. "+boost::lexical_cast<std::string>(pbc.cell[2][2])+" zlo zhi\n";
  //res += boost::lexical_cast<std::string>(pbc.cell[2][2])
  res += "\nMasses\n\n 1 55.85\n Atoms\n\n";
  return res;
};

void LAMMPSSimulator::lammps_path_write(std::string fn, double r) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<header();

  double ncom[]={0.,0.,0.};
  double c,nm=0.;

  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) {
    ncom[j] += pathway[3*i+j].deriv(1,r) / natoms;
  }

  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) {
    c = pathway[3*i+j].deriv(1,r)-ncom[j];
    nm += c * c;
  }
  nm = sqrt(nm);

  for(int i=0;i<natoms;i++){
    out<<i+1<<" 1 "; // TODO multispecies
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" "; // x y z
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" "; // path
    for(int j=0;j<3;j++) out<<(pathway[3*i+j].deriv(1,r)-ncom[j])/nm<<" ";
    for(int j=0;j<3;j++) out<<pathway[3*i+j].deriv(2,r)/nm/nm<<" ";
    out<<std::endl;
  }
  out.close();
};


void LAMMPSSimulator::lammps_dev_write(std::string fn, double r, double *dev, double *dev_sq) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<header();
  for(int i=0;i<natoms;i++){
    out<<i+1<<" 1 "; // TODO multispecies
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" "; // x y z
    for(int j=0;j<3;j++) out<<dev[3*i+j]<<" "; // <dev>
    for(int j=0;j<3;j++) out<<sqrt(dev_sq[3*i+j]-dev[3*i+j]*dev[3*i+j])<<" ";
    out<<std::endl;
  }
  out.close();
};
