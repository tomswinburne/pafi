#include "GeneralSimulator.hpp"

void GeneralSimulator::write(std::string fn, double r) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
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

  for (int i=0;i<natoms;i++) {
    out<<i<<" ";
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" "; // x y z
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" "; // path
    for(int j=0;j<3;j++) out<<(pathway[3*i+j].deriv(1,r)-ncom[j])/nm<<" ";
    for(int j=0;j<3;j++) out<<pathway[3*i+j].deriv(2,r)/nm/nm<<" ";
    out<<std::endl;
  }
  out.close();
};

void GeneralSimulator::write_dev(std::string fn, double r, double *dev, double *dev_sq) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<"# PAFI DUMP FILE. Reference path u(r) is a Nx3 vector.\n";
  out<<"# For i=0,1,2: u_i(r) , mean(x_i-u_i|r) , std(x_i-u_i|r) across ensemble\n";
  for(int i=0;i<natoms;i++) {
    out<<i+1<<" ";
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" ";
    for(int j=0;j<3;j++) out<<dev[3*i+j]<<" ";
    for(int j=0;j<3;j++) out<<dev_sq[3*i+j]-dev[3*i+j]*dev[3*i+j]<<" ";
    out<<std::endl;
  }
  out.close();
};


void GeneralSimulator::expansion(double T, double *newscale) {
  double coeff;

  newscale[0] = 1.0;
  coeff = std::stod(params->parameters["LinearThermalExpansionX"]);
  newscale[0] += coeff*T;
  coeff = std::stod(params->parameters["QuadraticThermalExpansionX"]);
  newscale[0] += coeff*T*T;

  newscale[1] = 1.0;
  coeff = std::stod(params->parameters["LinearThermalExpansionY"]);
  newscale[1] += coeff*T;
  coeff = std::stod(params->parameters["QuadraticThermalExpansionY"]);
  newscale[1] += coeff*T*T;

  newscale[2] = 1.0;
  coeff = std::stod(params->parameters["LinearThermalExpansionZ"]);
  newscale[2] += coeff*T;
  coeff = std::stod(params->parameters["QuadraticThermalExpansionZ"]);
  newscale[2] += coeff*T*T;
};

void GeneralSimulator::make_path(std::vector<std::string> knot_list, int real_coord) {
  pathway_r.clear();
  int nknots = knot_list.size();
  // no way around it- have to store all the knots
  std::vector<double> xs(nknots,0.), ys(nknots,0.), zs(nknots,0.);
  std::vector<double> r(nknots,0.), rr(nknots,0.);
  double *knots = new double[(const int)(3*natoms*nknots)];
  double *x = new double[(const int)(3*natoms)];
  spline xspl,yspl,zspl;
  double dx;

  // run through knots, and make spline
  load_config(knot_list[0],x);

  for (int i=0;i<3*natoms;i++) knots[i] = x[i];

  for (int knot=1; knot<nknots; knot++) {
    load_config(knot_list[knot],x);
    for(int i=0;i<3*natoms;i++) x[i]-=knots[i];
    pbc.wrap(x,3*natoms);
    for(int i=0;i<3*natoms;i++) \
      knots[i+knot*3*natoms] = x[i]+knots[i];
  }
  
  if(real_coord>0) {
    /*
      real_coord==1 : use the real space distance between knots
      real_coord==2 : use symmetrized distance
    */

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

    if(real_coord==1) {
      for(int knot=0;knot<nknots;knot++) pathway_r.push_back(r[knot]);
    } else {
      for(int knot=0;knot<nknots;knot++) {
        pathway_r.push_back(0.5*(r[knot] + 1.0 - rr[knot]));
        r[knot] = 0.5*(r[knot] + 1.0 - rr[knot]);
      }
    }    
  } else {
    for(int knot=0;knot<nknots;knot++) {
      r[knot] = 1.0*float(knot)/float(nknots-1);
      pathway_r.push_back(r[knot]);
    }
  }

  for(int i=0; i<natoms; i++) {
    for(int knot=0;knot<nknots;knot++) {
      xs[knot] = knots[3*natoms*knot + 3*i + 0];
      ys[knot] = knots[3*natoms*knot + 3*i + 1];
      zs[knot] = knots[3*natoms*knot + 3*i + 2];
    }


    xspl.set_points(r,xs,params->cubic_spline);
    pathway.push_back(xspl);

    yspl.set_points(r,ys,params->cubic_spline);
    pathway.push_back(yspl);

    zspl.set_points(r,zs,params->cubic_spline);
    pathway.push_back(zspl);

  }
  delete [] knots; // clear memory
};

double GeneralSimulator::path(int i, double r, int d, double s) {
  if(d==0) return pathway[i].deriv(d,r) * s;
  double dr = 1.0 / pathway_r.size();
  double val = pathway[i].deriv(0,r);
  if(d==1) return (pathway[i].deriv(0,r+dr)-val) * s/dr;
  else return (pathway[i].deriv(0,r+dr)+pathway[i].deriv(0,r-dr)-2*val) * s/dr/dr;
};

void GeneralSimulator::evaluate(std::vector<double> &results) {
  // Rescale, establish hp fix
  getEnergy();
};
