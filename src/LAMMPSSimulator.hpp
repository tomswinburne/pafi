#ifndef LSIM_H
#define LSIM_H

#include "GeneralSimulator.hpp"

class LAMMPSSimulator : public GeneralSimulator {

public:
  LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int rank) : GeneralSimulator (instance_comm, p, rank) {
    tag = rank;
    params = &p;
    scale = 1.;
    fixhp = false;
    temperature = 0.;
    position = 0.;
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

  virtual void load_config(std::string file_string,std::vector<double> &x) {
    lammps_command(lmp,(char *)"delete_atoms group all");
    std::string fn = "read_data ";
    fn += file_string;
    fn += " add merge";
    lammps_command(lmp,(char *)fn.c_str());
    lammps_gather_atoms(lmp,(char *) "x",1,3,&x[0]);
  };

  virtual void run_script(std::string sn){
    std::vector<std::string> strv = params->Script(sn);
    for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
  };

  virtual void run_commands(std::vector<std::string> strv) {
    for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
  };

  virtual void populate(double r) { // kind of LAMMPS dependent...

    std::vector<double> t(3*natoms,0.),s(3*natoms,0.);
    double ncom[]={0.,0.,0.};
    norm_mag=0.;

    // fill x and path
    for(int i=0;i<3*natoms;i++) t[i] = pathway[i](r) * scale;

    lammps_scatter_atoms(lmp,(char *)"x",1,3,&t[0]);

    lammps_command(lmp,(char *)"run 0"); // ESSENTIAL FOR REINDEXING!

    lammps_scatter_atoms(lmp,(char *)"path",1,3,&t[0]);
    // Fill normal==dpath and dnormal==ddpath. See paper for normalization
    for(int i=0;i<3*natoms;i++) {
      t[i] = pathway[i].deriv(1,r) * scale;
      ncom[i%3] += t[i]/(double)natoms;
    }
    for(int i=0;i<3*natoms;i++){
      t[i] -= ncom[i%3];
      norm_mag += t[i]*t[i];
    }

    norm_mag = sqrt(norm_mag);
    for(int i=0;i<3*natoms;i++) t[i] /= norm_mag;
    lammps_scatter_atoms(lmp,(char *)"norm",1,3,&t[0]);
    for(int i=0;i<3*natoms;i++) \
      t[i] = pathway[i].deriv(2,r) * scale / norm_mag / norm_mag;
    lammps_scatter_atoms(lmp,(char *)"dnorm",1,3,&t[0]);
    lammps_command(lmp,(char *)"run 0");

    /*
    for(int i=0;i<3*natoms;i++) t[i]=0.;
    lammps_gather_atoms(lmp,(char *)"norm",1,3,&t[0]);
    //lammps_command(lmp,(char *)"run 0");

    lammps_gather_atoms(lmp,(char *)"norm",1,3,&s[0]);

    double ttnn;
    for(int i=0;i<3*natoms;i++) {
      ttnn = (pathway[i].deriv(1,r)-ncom[i%3])/norm_mag;
      if(fabs(t[i]-ttnn)>0.001) std::cout<<"PROBLEM: "<<i<<" "<<t[i]<<" "<<ttnn<<" "<<ttnn-t[i]<<"\n";
    }

    std::ofstream out;
    std::string fn;
    fn = "PathNorm.dat";
    out.open(fn.c_str(),std::ofstream::out);
    for(int i=0;i<natoms;i++) {
      out<<i<<" ";
      for(int j=0;j<3;j++) out<<s[3*i+j]<<" ";
      for(int j=0;j<3;j++) out<<t[3*i+j]<<" ";
      for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" ";
      for(int j=0;j<3;j++) out<<(pathway[3*i+j].deriv(1,r)-ncom[j])/norm_mag<<" ";
      out<<"\n";
    }
    out.close();
    */


  };

  virtual void setup(double r, double T){
    //if(tag==0) std::cout<<"SETUP\n";
    position = r;
    temperature = T;

    double coeff,new_scale = 1.0;
    coeff = boost::lexical_cast<double>(params->parameters["Linear"]);
    new_scale += coeff*T;
    coeff = boost::lexical_cast<double>(params->parameters["Quadratic"]);
    new_scale += coeff*T*T;
    std::string cmd = "change_box all x scale ";
    cmd += boost::lexical_cast<std::string>(new_scale/scale);
    cmd += " y scale ";
    cmd += boost::lexical_cast<std::string>(new_scale/scale);
    cmd += " z scale ";
    cmd += boost::lexical_cast<std::string>(new_scale/scale);
    //cmd += " remap";
    lammps_command((void *)lmp, (char *)cmd.c_str());

    params->parameters["Temperature"] = boost::lexical_cast<std::string>(T);
    if (!fixhp) {
      //if(tag==0) std::cout<<"FHP"<<std::endl;
      cmd = "fix hp all hp %Temperature% 0.1 %RANDOM% overdamped 1 com 1\nrun 0";
      fixhp = true;
    }
    run_commands(params->Parse(cmd));
    scale = new_scale;
    populate(position);
    cmd = "run 0";
    run_commands(params->Parse(cmd));
    refE = getEnergy();
  };

  double thermalize() {
    double T;
    std::string cmd = "reset_timestep 0\n";
    cmd += "fix ape all ave/time 1 %TWindow% %TWindow% v_pe\nrun %ThermSteps%";
    run_commands(params->Parse(cmd));
    double * Tptr = (double *) lammps_extract_fix(lmp,(char *)"ape",0,0,0,0);
    T = (*Tptr-refE)/natoms/1.5/8.617e-5;
    cmd = "unfix ape";
    run_commands(params->Parse(cmd));
    lammps_free(Tptr);
    return T;
  };

  double sample(std::vector<double> &results) {
    double afd;
    std::string cmd = "reset_timestep 0\nrun %SampleSteps%";
    //cmd += "fix ape all ave/time 1 %TWindow% %TWindow% v_pe\n";
    //cmd += "fix af all ave/time 1 %SampleSteps% %SampleSteps% f_hp[1] f_hp[2] f_hp[3]\nrun %SampleSteps%";
    run_commands(params->Parse(cmd));
    //if (tag==0) std::cout<<"RAN HP\n";
    //double * Tptr = (double *) lammps_extract_fix(lmp,(char *)"ape",0,0,0,0);
    //double T = (*Tptr-refE)/natoms/1.5/8.617e-5;
    double * afptr = (double *) lammps_extract_fix(lmp,(char *)"hp",0,1,0,0);
    results.clear();
    afd = *afptr;
    results.push_back( afd * norm_mag );
    //cmd = "unfix af\nunfix ape";
    //cmd = "unfix ape";
    //run_commands(params->Parse(cmd));
    //lammps_free(afptr);
    //lammps_free(Tptr);
    return 0.;
    //return T;
  };

  double getEnergy(){
    lammps_command(lmp,(char *) "run 0");
    void * lmpE = lammps_extract_compute(lmp,(char *) "pe",0,0);
    if (lmpE == NULL){
      lammps_command(lmp,(char *) "compute pe all pe");
      lammps_command(lmp,(char *) "variable pe equal pe");
      lammps_command(lmp,(char *) "run 0");
      lmpE = lammps_extract_compute(lmp,(char *) "pe",0,0);
    }
    double baseE = *((double *) lmpE);
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
