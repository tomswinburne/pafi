#ifndef LSIM_H
#define LSIM_H

#include "GeneralSimulator.hpp"

class LAMMPSSimulator : public GeneralSimulator {

public:
  LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int rank) : GeneralSimulator (instance_comm, p, rank) {
    tag = rank;
    params = &p;

    // set up LAMMPS
    char str1[32];
    char **lmparg = new char*[5];
    lmparg[0] = NULL; // required placeholder for program name
    lmparg[1] = (char *) "-screen";
    lmparg[2] = (char *) "none";
    lmparg[3] = (char *) "-log";
    lmparg[4] = (char *) "none";
    sprintf(str1,"log.lammps.%d.seule",rank);
    if (rank==0) lmparg[4] = str1;

    lammps_open(5,lmparg,instance_comm,(void **) &lmp);
    run_script("Input");
    natoms = *((int *) lammps_extract_global(lmp,(char *) "natoms"));

    id = std::vector<int>(natoms,0);
    lammps_gather_atoms(lmp,(char *) "id",0,1,&id[0]);

    // get cell info
    pbc.load(lmp);

    // Get type / image info
    species = std::vector<int>(natoms,1);
    lammps_gather_atoms(lmp,(char *) "type",0,1,&species[0]);

    image = std::vector<int>(natoms,1);
    lammps_gather_atoms(lmp,(char *) "image",0,1,&image[0]);

    // prepare for fix?

  };

  /*
    Load xyz configuration from data file and put in vector
  */
  void load_config(std::string file_string,std::vector<double> &x) {
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
  void run_script(std::string sn){
    std::vector<std::string> strv = params->Script(sn);
    for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
  };

  /*
    Parse and run script from string with linebreaks
  */
  void run_commands(std::vector<std::string> strv) {
    for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
  };

  /*
    Fill configuration, path, tangent and tangent gradient. Return tangent norm
  */
  void populate(double r, double scale, double &norm_mag) {
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
  void rescale_cell(double scale) {
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
  void sample(double r, double T, std::vector<double> &results, std::vector<double> &deviation) {
    std::string cmd;
    double norm_mag, scale, sampleT;
    double *lmp_ptr;
    double **lmp_dev_ptr;

    scale = expansion(T);
    rescale_cell(scale);
    populate(r,scale,norm_mag);
    //if(tag==0) std::cout<<"POPULATED\n";

    deviation.clear();
    //lmp_dev_ptr = (double **) lammps_extract_fix(lmp,(char *)"ax",1,2,0,0);
    for(int i=0;i<3*natoms;i++) deviation.push_back(0.);
    lammps_gather_atoms(lmp,(char *) "x",1,3,&deviation[0]);
    for(int i=0;i<3*natoms;i++) deviation[i] -= pathway[i](r);
    double dm=0.;
    pbc.wrap(deviation);
    for(int i=0;i<3*natoms;i++) dm += deviation[i] * deviation[i];
    if(tag==0) std::cout<<dm<<std::endl;

    params->parameters["Temperature"] = boost::lexical_cast<std::string>(T);
    cmd = "fix hp all hp %Temperature% 0.01 %RANDOM% overdamped 1 com 1\nrun 0";
    run_commands(params->Parse(cmd));

    refE = getEnergy();
    cmd = "reset_timestep 0\n";
    cmd += "fix ae all ave/time 1 %TWindow% %TWindow% v_pe\n";
    cmd += "run %ThermSteps%";
    run_commands(params->Parse(cmd));

    results.clear();
    lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
    sampleT = (*lmp_ptr-refE)/natoms/1.5/8.617e-5;
    cmd = "unfix ae\ncompute xyzu all property/atom xu yu zu\nrun 0";
    run_commands(params->Parse(cmd));
    results.push_back(sampleT);

    cmd = "reset_timestep 0\n";
    cmd += "fix ae all ave/time 1 %TWindow% %TWindow% v_pe\n";
    //cmd += "fix ad all ave/deviation 1 %TWindow% %TWindow%\n";
    cmd += "fix axu all ave/atom 1 %TWindow% %TWindow% c_xyzu[1] c_xyzu[2] c_xyzu[2]\n";
    cmd += "fix af all ave/time 1 %SampleSteps% %SampleSteps% f_hp[1] f_hp[2] f_hp[3]\nrun %SampleSteps%";
    run_commands(params->Parse(cmd));

    lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"ae",0,0,0,0);
    sampleT = (*lmp_ptr-refE)/natoms/1.5/8.617e-5;
    results.push_back(sampleT);

    lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,0,0);
    results.push_back(*lmp_ptr * norm_mag);

    lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,1,0);
    results.push_back(*lmp_ptr * norm_mag * norm_mag - results[2] * results[2]);

    lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,2,0);
    results.push_back(*lmp_ptr);

    lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"af",0,1,3,0);
    results.push_back(*lmp_ptr);

    //lammps_free(lmp_ptr);

    // TODO dump path efficiently....
    deviation.clear();
    for(int i=0;i<3*natoms;i++) deviation.push_back(0.);
    //lammps_gather_atoms(lmp,(char *) "x",1,3,&deviation[0]);
    //lmp_dev_ptr = (double **) lammps_extract_fix(lmp,(char *)"axu",1,2,0,0);

    for(int i=0;i<natoms;i++) for(int j=0;j<3;j++) {
      lmp_ptr = (double *) lammps_extract_fix(lmp,(char *)"axu",1,2,i,j);
      deviation[3*i+j] = *lmp_ptr-pathway[3*i+j](r);
    }
    dm=0.;
    //pbc.wrap(deviation);
    for(int i=0;i<3*natoms;i++) dm += deviation[i] * deviation[i];
    if(tag==0) std::cout<<dm<<std::endl;

    cmd = "unfix ae\nunfix af\nunfix axu\nunfix hp\nuncompute xyzu";
    run_commands(params->Parse(cmd));

    scale = 1. / scale;
    populate(r,scale,norm_mag);
    rescale_cell(scale);
  };

  double getEnergy() {
    lammps_command(lmp,(char *) "run 0");
    double * lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
    if (lmpE == NULL) {
      std::cout<<"HERE"<<std::endl;
      lammps_command(lmp,(char *) "compute pe all pe");
      lammps_command(lmp,(char *) "variable pe equal pe");
      lammps_command(lmp,(char *) "run 0");
      lmpE = (double *) lammps_extract_compute(lmp,(char *) "pe",0,0);
    }
    double baseE = *lmpE;
    //std::cout<<"HERE"<<std::endl;
    return baseE;
  };

  void close() {
    lammps_close(lmp);
  };

  // LAMMPS specific
  std::vector<int> species,q,image,id;
private:
  void *lmp;
};

#endif
