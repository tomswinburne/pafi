#ifndef LSYS_H
#define LSYS_H

#include "lammps/lammps.h"
#include "lammps/input.h"
#include "lammps/atom.h"
#include "lammps/library.h"

#include "spline.hpp"
#include "pbc.hpp"
#include "parser.hpp"

// To be templated.....

class LAMMPSSimulator {

public:
  LAMMPSSimulator (MPI_Comm &instance_comm, Parser &p, int rank) {
    tag = rank;
    params = &p;
    scale = 1.;
    norm_mag=1.;
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

    // get cell info
    pbc.load(lmp);

    // Get type / image info
    species = std::vector<int>(natoms,1);
    lammps_gather_atoms(lmp,(char *) "type",0,1,&species[0]);

    image = std::vector<int>(natoms,1);
    lammps_gather_atoms(lmp,(char *) "image",0,1,&image[0]);

    // prepare for fix?

  };

  void load_config(std::string file_string,std::vector<double> &x) {
    lammps_command(lmp,(char *)"delete_atoms group all");
    std::string fn = "read_data ";
    fn += file_string;
    fn += " add merge";
    lammps_command(lmp,(char *)fn.c_str());
    lammps_gather_atoms(lmp,(char *) "x",1,3,&x[0]);
  };

  void fill_vector(std::string vn, std::vector<double> &v) {
    lammps_scatter_atoms(lmp,(char *) vn.c_str(),1,3,&v[0]);
  };
  void fill_vector(std::string vn, std::vector<int> &v) {
    lammps_scatter_atoms(lmp,(char *) vn.c_str(),0,1,&v[0]);
  };
  /*double temperature() {
    return (double)*((double *) lammps_extract_compute(lmp,(char *) "thermo_temp",0,0));
  };*/

  void run_script(std::string sn){
    std::vector<std::string> strv = params->Script(sn);
    for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
  };

  void run_commands(std::vector<std::string> strv) {
    for(auto s:strv) lammps_command((void *)lmp,(char *)s.c_str());
  };


  void setup(double r, double T){
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
      std::cout<<"FHP\n";
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
    std::string cmd = "reset_timestep 0\n";
    cmd += "fix ape all ave/time 1 %TWindow% %TWindow% v_pe\n";
    cmd += "fix af all ave/time 1 %SampleSteps% %SampleSteps% f_hp[1] f_hp[2] f_hp[3]\nrun %SampleSteps%";
    run_commands(params->Parse(cmd));
    //if (tag==0) std::cout<<"RAN HP\n";
    //double * Tptr = (double *) lammps_extract_fix(lmp,(char *)"ape",0,0,0,0);
    //double T = (*Tptr-refE)/natoms/1.5/8.617e-5;
    //double * afptr = (double *) lammps_extract_fix(lmp,(char *)"hp",0,1,0,0);
    //results.clear();
    //afd = *afptr;
    //results.push_back( afd * norm_mag );
    cmd = "unfix af\nunfix ape";
    //cmd = "unfix ape";
    run_commands(params->Parse(cmd));
    //lammps_free(afptr);
    //lammps_free(Tptr);
    //return T;
    return 0.;
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

  // EVERYTHING BELOW SHOULD BE ENGINE INDEPENDENT

  void make_path(std::vector<std::string> knot_list) {
    // run through knots, and make spline
    int nknots = knot_list.size();
    std::vector<double> x(3*natoms);
    load_config(knot_list[0],x);
    // no way around it- have to store all the knots
    std::vector<double> xs(nknots,0.), ys(nknots,0.), zs(nknots,0.);
    std::vector<double> r(nknots,0.), rr(nknots,0.), knots(3*natoms*nknots,0.);

    for(int i=0;i<3*natoms;i++) knots[i] = x[i];
    for(int knot=1;knot<nknots;knot++) {
      load_config(knot_list[knot],x);
      for(int i=0;i<3*natoms;i++) x[i]-=knots[i];
      pbc.wrap(x);
      for(int i=0;i<3*natoms;i++) {
        knots[i+knot*3*natoms] = x[i]+knots[i];
      }
    }
    double dx;
    for(int knot=0;knot<nknots;knot++) {
      r[knot] = 0.;
      rr[knot] = 0.;
      for(int i=0;i<3*natoms;i++) {
        dx = knots[i+knot*3*natoms]-knots[i];
        r[knot] += dx*dx;
        dx = knots[i+knot*3*natoms]-knots[i+(nknots-1)*3*natoms];
        rr[knot] += dx*dx;
      }
    }
    for(int knot=0;knot<nknots-1;knot++) r[knot] = sqrt(r[knot]/r[nknots-1]);
    for(int knot=1;knot<nknots;knot++) rr[knot] = sqrt(rr[knot]/rr[0]);
    rr[0] = 1.0;
    r[nknots-1] = 1.0;
    for(int knot=0;knot<nknots;knot++) r[knot] = 0.5*r[knot]+0.5-0.5*rr[knot];
    for(int i=0; i<natoms; i++) {
      for(int knot=0;knot<nknots;knot++) {
        xs[knot] = knots[3*natoms*knot + 3*i+0];
        ys[knot] = knots[3*natoms*knot + 3*i+1];
        zs[knot] = knots[3*natoms*knot + 3*i+2];
      }
      spline::spline xspl,yspl,zspl;

      xspl.set_points(r,xs);
      pathway.push_back(xspl);

      yspl.set_points(r,ys);
      pathway.push_back(yspl);

      zspl.set_points(r,zs);
      pathway.push_back(zspl);
    }
    // delete knot array
    knots.clear();
    knots.shrink_to_fit();
    //load_config(knot_list[0],x);
  };

  void populate(double r) {
    std::vector<double> t(3*natoms,0.);
    double ncom[]={0.,0.,0.};
    norm_mag=0.;

    // fill x and path
    for(int i=0;i<3*natoms;i++) t[i] = pathway[i](r) * scale;

    fill_vector("x",t);
    fill_vector("path",t);

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
    fill_vector("norm",t);

    for(int i=0;i<3*natoms;i++) \
      t[i] = pathway[i].deriv(2,r) * scale / norm_mag / norm_mag;
    fill_vector("dnorm",t);
  };

  void evaluate(std::vector<double> &results) {
    // Rescale, establish hp fix
    std::cout<<"ENERGY: ";
    getEnergy();
  };
  double scale, position, temperature, refE, norm_mag;
  int natoms, tag;
  bool fixhp;
  MinImage pbc;
  Parser *params;
  std::vector<std::vector<double>> knots;
  std::vector<int> species,q,image;
  std::vector<spline::spline> pathway;
private:
  void *lmp;
};

#endif
