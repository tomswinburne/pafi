#include "GeneralSimulator.hpp"

/*void GeneralSimulator::write(double r,std::string fn) {
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

  for(int i=0;i<natoms;i++){
    out<<i+1<<" 1 "; // TODO multispecies
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" "; // x y z
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" "; // path
    for(int j=0;j<3;j++) out<<(pathway[3*i+j].deriv(1,r)-ncom[j])/nm<<" ";
    for(int j=0;j<3;j++) out<<pathway[3*i+j].deriv(2,r)/nm/nm<<" ";
    out<<std::endl;
  }
  out.close();
};*/


double GeneralSimulator::expansion(double T) {
  double coeff,new_scale = 1.0;
  coeff = boost::lexical_cast<double>(params->parameters["Linear"]);
  new_scale += coeff*T;
  coeff = boost::lexical_cast<double>(params->parameters["Quadratic"]);
  new_scale += coeff*T*T;
  return new_scale;
};

void GeneralSimulator::make_path(std::vector<std::string> knot_list) {
  int nknots = knot_list.size();
  // no way around it- have to store all the knots
  std::vector<double> x(3*natoms,0.);
  std::vector<double> xs(nknots,0.), ys(nknots,0.), zs(nknots,0.);
  std::vector<double> r(nknots,0.), rr(nknots,0.);
  double *knots = new double[(const int)(3*natoms*nknots)];
  spline xspl,yspl,zspl;
  double dx;

  // run through knots, and make spline
  load_config(knot_list[0],x);

  for (int i=0;i<3*natoms;i++) knots[i] = x[i];

  for (int knot=1; knot<nknots; knot++) {
    load_config(knot_list[knot],x);
    for(int i=0;i<3*natoms;i++) x[i]-=knots[i];
    pbc.wrap(x);
    for(int i=0;i<3*natoms;i++) \
      knots[i+knot*3*natoms] = x[i]+knots[i];
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

  for(int knot=0;knot<nknots;knot++) r[knot] = 0.5*(r[knot] + 1.0 - rr[knot]);

  for(int i=0; i<natoms; i++) {
    for(int knot=0;knot<nknots;knot++) {
      xs[knot] = knots[3*natoms*knot + 3*i + 0];
      ys[knot] = knots[3*natoms*knot + 3*i + 1];
      zs[knot] = knots[3*natoms*knot + 3*i + 2];
    }
    
    xspl.set_points(r,xs);
    pathway.push_back(xspl);

    yspl.set_points(r,ys);
    pathway.push_back(yspl);

    zspl.set_points(r,zs);
    pathway.push_back(zspl);

  }
  delete [] knots; // clear memory
};

void GeneralSimulator::evaluate(std::vector<double> &results) {
  // Rescale, establish hp fix
  getEnergy();
};
