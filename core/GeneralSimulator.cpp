#include "GeneralSimulator.hpp"

GeneralSimulator::GeneralSimulator(MPI_Comm &instance_comm, Parser &p, int t) {
  tag = t;
  comm = &instance_comm;
  MPI_Comm_rank(*comm,&local_rank);
  MPI_Comm_size(*comm,&local_size);
  parser = &p;
  error_count = 0;
  scale[0]=1.0; scale[1]=1.0; scale[2]=1.0;
  last_error_message="";
  out_width = 16;
  // to be overwritten
  natoms = 0;
  nlocal = 0;
  offset = 0;

};

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

void GeneralSimulator::write_dev(std::string fn, double r, double *dev) {
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<"# PAFI DUMP FILE. Reference path u(r) is a Nx3 vector.\n";
  out<<"# For i=0,1,2: u_i(r) , mean(x_i-u_i|r) across valid ensemble\n";
  for(int i=0;i<natoms;i++) {
    out<<i+1<<" ";
    for(int j=0;j<3;j++) out<<pathway[3*i+j](r)<<" ";
    for(int j=0;j<3;j++) out<<dev[3*i+j]<<" ";
    out<<std::endl;
  }
  out.close();
};


void GeneralSimulator::expansion(double T, double *newscale) {
  double coeff;

  newscale[0] = 1.0;
  coeff = std::stod(parser->parameters["LinearThermalExpansionX"]);
  newscale[0] += coeff*T;
  coeff = std::stod(parser->parameters["QuadraticThermalExpansionX"]);
  newscale[0] += coeff*T*T;

  newscale[1] = 1.0;
  coeff = std::stod(parser->parameters["LinearThermalExpansionY"]);
  newscale[1] += coeff*T;
  coeff = std::stod(parser->parameters["QuadraticThermalExpansionY"]);
  newscale[1] += coeff*T*T;

  newscale[2] = 1.0;
  coeff = std::stod(parser->parameters["LinearThermalExpansionZ"]);
  newscale[2] += coeff*T;
  coeff = std::stod(parser->parameters["QuadraticThermalExpansionZ"]);
  newscale[2] += coeff*T*T;
  //std::cout<<scale[0]<<" "<<scale[1]<<" "<<scale[2]<<std::endl;
};

void GeneralSimulator::make_path(std::vector<std::string> knot_list) {
  /*
    TODO: parallel i/o and splining
    Only really a problem with memory limitations, say 2GB / core.
    This implies 300M coordinates at double precision == 1M atoms, 100 planes.
    Typical large-scale use - 150k atoms, 20 planes.
    So memory-heavy but nothing problematic so far, leaving for future
  */
  pathway_r.clear();
  pathway.clear();

  if(nlocal==0) {
    if(local_rank==0)
      std::cout<<"GeneralSimulator::make_path : Not initialized!"<<std::endl;
    return ;
  }

  int nknots = knot_list.size();
  double dx;
  std::vector<double> xs(nknots,0.), r(nknots,0.), rr(nknots,0.);

  // big allocations, but have to do it at some point... is deleted before simulations
  double *knots = new double[nlocal*nknots];
  pathway.assign(nlocal,*(new spline));

  // run through knots, and make spline
  // in current implementation, nlocal = 3*natoms, offset = 0;
  load_config(knot_list[0],x);

  for (int i=0;i<nlocal;i++) knots[i] = x[i+offset];

  for (int knot=1; knot<nknots; knot++) {
    load_config(knot_list[knot],x);
    for(int i=0;i<nlocal;i++) x[i+offset]-=knots[i];
    pbc.wrap(x+offset,nlocal);
    for(int i=0;i<nlocal;i++) \
      knots[i+knot*nlocal] = x[i+offset]+knots[i];
  }

  for(int knot=0;knot<nknots;knot++) {
    r[knot] = 0.;
    rr[knot] = 0.;
    for(int i=0;i<nlocal;i++) {
      dx = knots[i+knot*nlocal]-knots[i];
      r[knot] += dx*dx;
      dx = knots[i+knot*nlocal]-knots[i+(nknots-1)*nlocal];
      rr[knot] += dx*dx;
    }
  }

  for(int knot=0;knot<nknots-1;knot++) r[knot] = sqrt(r[knot]/r[nknots-1]);
  for(int knot=1;knot<nknots;knot++) rr[knot] = sqrt(rr[knot]/rr[0]);
  rr[0] = 1.0;
  r[nknots-1] = 1.0;

  for(int knot=0;knot<nknots;knot++) {
    pathway_r.push_back(0.5*(r[knot] + 1.0 - rr[knot]));
    r[knot] = 0.5*(r[knot] + 1.0 - rr[knot]);
  }


  for(int i=0; i<nlocal; i++) {
    for(int knot=0;knot<nknots;knot++) xs[knot] = knots[nlocal*knot + i];
    pathway[i].set_points(r,xs,parser->spline_path);
  }
  delete [] knots; // clear memory

  sample_r.clear();
  double dr = 0.1;
  if (parser->nPlanes>1)
    dr = (parser->stopr-parser->startr)/(double)(parser->nPlanes-1);

  if(parser->spline_path and not parser->match_planes) {
    for (double r = parser->startr; r <= parser->stopr+0.5*dr; r += dr )
      sample_r.push_back(r);
  } else {
    for(auto r: pathway_r) if(r>=0.0 && r<=1.0) sample_r.push_back(r);
  }

};

double GeneralSimulator::path(int i, double r, int d, double s) {

  if(parser->spline_path or d==0) return pathway[i].deriv(d,r) * s;
  double dr = 1.0 / pathway_r.size();
  double val = pathway[i].deriv(0,r);
  if(d==1) return (pathway[i].deriv(0,r+dr)-val) * s/dr;
  else return (pathway[i].deriv(0,r+dr)+pathway[i].deriv(0,r-dr)-2*val) * s/dr/dr;
};

void GeneralSimulator::evaluate(std::vector<double> &results) {
  // Rescale, establish hp fix
  getEnergy();
};
