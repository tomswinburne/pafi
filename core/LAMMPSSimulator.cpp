#include "LAMMPSSimulator.hpp"

LAMMPSSimulator::LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int rank) {
  tag = rank;
  params = &p;
  scale = new double[3];
  // set up LAMMPS
  char str1[32];
  char **lmparg = new char*[5];
  lmparg[0] = NULL; // required placeholder for program name
  lmparg[1] = (char *) "-screen";
  lmparg[2] = (char *) "none";
  lmparg[3] = (char *) "-log";
  if(params->loglammps) {
    sprintf(str1,"log.lammps.%d",rank);
    lmparg[4] = str1;
  } else lmparg[4] = (char *) "none";

  lammps_open(5,lmparg,instance_comm,(void **) &lmp);
  run_script("Input");
  natoms = *((int *) lammps_extract_global(lmp,(char *) "natoms"));

  has_pafi = (bool)lammps_config_has_package((char *)"USER-MISC");

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

  made_fix=false;
  made_compute=false;

  expansion(0.0);
  rescale_cell();

  // prepare for fix?
  //lammps_command(lmp,"fix pafipath all property/atom d_nx d_ny d_nz d_dnx d_dny d_dnz d_ddnx d_ddny d_ddnz");
  //lammps_command(lmp,"run 0");
  //lammps_command(lmp,"compute pafipath all property/atom d_nx d_ny d_nz d_dnx d_dny d_dnz d_ddnx d_ddny d_ddnz");
};

/*
  Load xyz configuration from data file and put in vector
*/
void LAMMPSSimulator::load_config(std::string file_string,std::vector<double> &x) {
  made_fix=false;
  made_compute=false;
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
  for(auto s:strv) {
    #ifdef VERBOSE
    std::cout<<s<<std::endl;
    #endif
    lammps_command((void *)lmp,(char *)s.c_str());
  }
};

/*
  Parse and run script from string with linebreaks
*/
void LAMMPSSimulator::run_commands(std::vector<std::string> strv) {
  for(auto s:strv) {
    #ifdef VERBOSE
    std::cout<<s<<std::endl;
    #endif
    lammps_command((void *)lmp,(char *)s.c_str());
  }
};

void LAMMPSSimulator::run_commands(std::string strv) {
  run_commands(params->Parse(strv));
};


/*
  Fill configuration, path, tangent and tangent gradient. Return tangent norm
*/
void LAMMPSSimulator::populate(double r, double &norm_mag) {
  //std::cout<<scale[0]<<" "<<scale[1]<<" "<<scale[2]<<std::endl;
  std::vector<double> t(3*natoms,0.),lt(natoms,0.0);

  double ncom[]={0.,0.,0.};
  norm_mag=0.;
  // check for __pafipath fix and compute
  if (!made_fix) {
    #ifdef VERBOSE
    std::cout<<"making __pafipath fix"<<std::endl;
    #endif
    lammps_command(lmp,(char *)"fix __pafipath all property/atom d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz");
    lammps_command(lmp,(char *)"run 0");
    made_fix=true;
  }


  //int dtype;
  //int indexval = lmp->atom->find_custom("d_nx",dtype);
  //std::cout<<"INDEX: "<<indexval<<" "<<dtype<<std::endl;

  //int icompute = lmp->modify->find_compute("__pafipath");
  //if(icompute<0) {
  //  lammps_command(lmp,"compute pafipath all property/atom d_nx d_ny d_nz d_dnx d_dny d_dnz d_ddnx d_ddny d_ddnz");
  //  lammps_command(lmp,"run 0"); // not needed I think...
  //}

  std::string xyz[3]; xyz[0]="x"; xyz[1]="y"; xyz[2]="z";
  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) t[3*i+j] = pathway[3*i+j].deriv(0,r) * scale[j];
  lammps_scatter_atoms(lmp,(char *)"x",1,3,&t[0]);

  for(int j=0;j<3;j++) {
    for(int i=0;i<natoms;i++)  lt[i] = pathway[3*i+j].deriv(0,r) * scale[j];
    #ifdef VERBOSE
      std::cout<<"Scattering d_u"+xyz[j]<<std::endl;
    #endif
    lammps_scatter_atoms(lmp,(char *)("u"+xyz[j]).c_str(),1,1,&lt[0]);
  }

  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) t[3*i+j] = pathway[3*i+j].deriv(1,r) * scale[j];

  // Center of mass projection and tangent length
  for(int i=0;i<3*natoms;i++) ncom[i%3] += t[i]/(double)natoms;
  for(int i=0;i<3*natoms;i++) t[i] -= ncom[i%3];
  for(int i=0;i<3*natoms;i++) norm_mag += t[i]*t[i];
  norm_mag = sqrt(norm_mag);
  for(int i=0;i<3*natoms;i++) t[i] /= norm_mag;
  for(int j=0;j<3;j++) {
    for(int i=0;i<natoms;i++)  lt[i] = t[3*i+j];
    lammps_scatter_atoms(lmp,(char *)("n"+xyz[j]).c_str(),1,1,&lt[0]);
    for(int i=0;i<natoms;i++)  lt[i] = pathway[3*i+j].deriv(2,r) * scale[j] / norm_mag / norm_mag;
    lammps_scatter_atoms(lmp,(char *)("dn"+xyz[j]).c_str(),1,1,&lt[0]);
  }

  lammps_command(lmp,(char *)"run 0");

  if(!made_compute) {
    #ifdef VERBOSE
    std::cout<<"making __pafipath compute"<<std::endl;
    #endif
    lammps_command(lmp,(char *)"compute __pafipath all property/atom d_ux d_uy d_uz d_nx d_ny d_nz d_dnx d_dny d_dnz");
    made_compute=true;
    lammps_command(lmp,(char *)"run 0");
    //lammps_command(lmp,(char *)"run 0");
  }
};

/*
  Rescale simulation cell
*/
void LAMMPSSimulator::rescale_cell() {
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
  sample temperature <f>, <f^2>, <psi> and <x-u>.n, and max_jump
*/
void LAMMPSSimulator::sample(double r, double T, double *results, double *dev) {
  std::string cmd;
  double norm_mag, sampleT, dm;
  double *lmp_ptr;
  double **lmp_dev_ptr;
  for(int j=0;j<3;j++) scale[j] = 1.0;
  populate(r,norm_mag);

  // Stress Fixes
  run_script("PreRun");
  lammps_command(lmp,(char *)"run 0");

  expansion(T);
  rescale_cell();
  populate(r,norm_mag);


  params->parameters["Temperature"] = std::to_string(T);
  cmd = "fix hp all pafi __pafipath "+params->parameters["Temperature"]+" ";
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
  sampleT = (*lmp_ptr-refE)/natoms/1.5/BOLTZ;
  lammps_free(lmp_ptr);
  results[0] = sampleT;

  cmd = "unfix ae\nrun 0";
  run_commands(cmd);
  //results[0] = sampleT;

  // TODO: groupname for ave/atom
  std::string SampleSteps = params->parameters["SampleSteps"];
  cmd = "reset_timestep 0\n";
  cmd += "fix ae all ave/time 1 "+SampleSteps+" "+SampleSteps+" v_pe\n";
  cmd += "fix ap all ave/atom 1 "+SampleSteps+" "+SampleSteps+" x y z\n";
  cmd += "fix af all ave/time 1 "+SampleSteps+" "+SampleSteps+" f_hp[1] f_hp[2] f_hp[3] f_hp[4]\n";
  cmd += "run "+SampleSteps;
  run_commands(cmd);

  lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
  sampleT = (*lmp_ptr-refE)/natoms/1.5/BOLTZ;
  lammps_free(lmp_ptr);
  results[1] = sampleT;

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


  cmd = "min_style fire\n minimize 0 0 ";
  cmd += params->parameters["ThermSteps"]+" "+params->parameters["ThermSteps"];
  run_commands(cmd);
  lammps_gather_atoms(lmp,(char *)"x",1,3,dev);
  for(int i=0;i<3*natoms;i++) dev[i] -= pathway[i].deriv(0,r)*scale[i%3];
  pbc.wrap(dev,3*natoms);
  double a_disp=0.0,max_disp = 0.0;
  for(int i=0;i<natoms;i++) {
    a_disp = 0.0;
    for(int j=0;j<3;j++) a_disp += dev[3*i+j]*dev[3*i+j];
    max_disp = std::max(a_disp,max_disp);
    for(int j=0;j<3;j++) dev[3*i+j] = 0.0;
  }
  
  lammps_gather_fix(lmp,(char *)"ap",1,3,dev);
  for(int i=0;i<3*natoms;i++) dev[i] -= pathway[i].deriv(0,r)*scale[i%3];
  pbc.wrap(dev,3*natoms);



  dm=0.; for(int i=0;i<3*natoms;i++) dm += dev[i] * dev[i];
  results[6] = dm;

  results[7] = max_disp;

  cmd = "unfix ae\nunfix af\nunfix ap\nunfix hp";
  run_commands(cmd);

  // Stress Fixes
  run_script("PostRun");
  cmd = "run 0"; run_commands(cmd);

  // rescale back
  for(int j=0;j<3;j++) scale[j] = 1.0/scale[j];
  populate(r,norm_mag);
  rescale_cell();
  for(int j=0;j<3;j++) scale[j] = 1.0;
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

std::string LAMMPSSimulator::header(double mass=55.85) {
  std::string res="LAMMPS dump file\n\n";
  res += std::to_string(natoms)+" atoms\n";
  res += "1 atom types\n\n"; // TODO species
  res += "0. "+std::to_string(pbc.cell[0][0]*scale[0])+" xlo xhi\n";
  res += "0. "+std::to_string(pbc.cell[1][1]*scale[1])+" ylo yhi\n";
  res += "0. "+std::to_string(pbc.cell[2][2]*scale[2])+" zlo zhi\n";
  res += std::to_string(pbc.cell[0][1]*scale[1])+" ";
  res += std::to_string(pbc.cell[0][2]*scale[2])+" ";
  res += std::to_string(pbc.cell[1][2]*scale[2])+" xy xz yz\n";
  //res += std::to_string(pbc.cell[2][2])
  res += "\nMasses\n\n 1 "+std::to_string(mass)+"\n Atoms\n\n";
  return res;
};

void LAMMPSSimulator::lammps_dump_path(std::string fn, double r) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<header();

  double ncom[]={0.,0.,0.};
  double c,nm=0.;

  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) {
    ncom[j] += pathway[3*i+j].deriv(1,r)*scale[j] / natoms;
  }

  for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) {
    c = pathway[3*i+j].deriv(1,r)*scale[j]-ncom[j];
    nm += c * c;
  }
  nm = sqrt(nm);

  for(int i=0;i<natoms;i++){
    out<<i+1<<" 1 "; // TODO multispecies
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)*scale[j]<<" "; // x y z
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)*scale[j]<<" "; // path
    for(int j=0;j<3;j++) out<<(pathway[3*i+j].deriv(1,r)*scale[j]-ncom[j])/nm<<" ";
    for(int j=0;j<3;j++) out<<pathway[3*i+j].deriv(2,r)*scale[j]/nm/nm<<" ";
    out<<std::endl;
  }
  out.close();
};

void LAMMPSSimulator::lammps_write_dev(std::string fn, double r, double *dev, double *dev_sq) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<header();
  for(int i=0;i<natoms;i++){
    out<<i+1<<" 1 "; // TODO multispecies
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)*scale[j]<<" "; // x y z
    for(int j=0;j<3;j++) out<<dev[3*i+j]<<" "; // <dev>
    for(int j=0;j<3;j++) out<<sqrt(dev_sq[3*i+j]-dev[3*i+j]*dev[3*i+j])<<" ";
    out<<std::endl;
  }
  out.close();
};
