#ifndef SIM_H
#define SIM_H

class GeneralSimulator {

public:
  GeneralSimulator (MPI_Comm &instance_comm, Parser &p, int rank) {
    tag = rank;
    params = &p;
    scale = 1.;
    fixhp = false;
    // set up SIMULATOR, pbc
  };

  virtual void load_config(std::string file_string,std::vector<double> &x){};

  virtual void run_script(std::string sn){};

  virtual void run_commands(std::vector<std::string> strv){};

  virtual void setup(double r, double T){};

  virtual double thermalize(){};

  virtual double sample(std::vector<double> &results){};

  virtual double getEnergy(){};

  virtual void close(){};

  virtual void populate(double r){};

  // LAMMPS INDEPENDENT
  virtual void make_path(std::vector<std::string> knot_list) {
    int nknots = knot_list.size();
    // no way around it- have to store all the knots
    std::vector<double> x(3*natoms);
    std::vector<double> xs(nknots,0.), ys(nknots,0.), zs(nknots,0.);
    std::vector<double> r(nknots,0.), rr(nknots,0.);
    std::vector<double> knots(3*natoms*nknots,0.);
    spline::spline xspl,yspl,zspl;
    double dx;

    // run through knots, and make spline
    load_config(knot_list[0],x);

    for(int i=0;i<3*natoms;i++) knots[i] = x[i];

    for(int knot=1;knot<nknots;knot++) {
      load_config(knot_list[knot],x);
      for(int i=0;i<3*natoms;i++) x[i]-=knots[i];
      pbc.wrap(x);
      for(int i=0;i<3*natoms;i++) {
        knots[i+knot*3*natoms] = x[i]+knots[i];
      }
    }

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

      xspl.set_points(r,xs);
      pathway.push_back(xspl);

      yspl.set_points(r,ys);
      pathway.push_back(yspl);

      zspl.set_points(r,zs);
      pathway.push_back(zspl);
    }
    knots.clear(); // delete knot array
  };

  virtual void evaluate(std::vector<double> &results) {
    // Rescale, establish hp fix
    getEnergy();
  };

  double scale, position, temperature, refE, norm_mag;
  int natoms, tag, nknots;
  bool fixhp;
  MinImage pbc;
  Parser *params;
  std::vector<spline::spline> pathway;
private:
  /* nothing */
};

#endif
