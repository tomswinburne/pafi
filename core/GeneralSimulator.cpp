#include "GeneralSimulator.hpp"

GeneralSimulator::GeneralSimulator(MPI_Comm &instance_comm, Parser &p, int rank) {
  position=0.0;
  tag = rank;
  comm = &instance_comm;
  params = &p;
  MPI_Comm_rank(*comm,&local_rank);
  MPI_Comm_size(*comm,&local_size);
};

void GeneralSimulator::write_dev(std::string fn, double r, double *dev, double *dev_sq) {
  int i,j;
  std::ofstream out;
  out.open(fn.c_str(),std::ofstream::out);
  out<<"# PAFI DUMP FILE. Reference path u(r) is a Nx3 vector.\n";
  out<<"# For i=0,1,2: u_i(r) , < x_i-u_i | r > , <(x_i-u_i)^2 | r >)\n";
  fill_path(r,0,glo_data);
  for(i=0;i<natoms;i++) {
    out<<i<<" ";
    for(j=0;j<3;j++) out<<glo_data[3*i+j]<<" ";
    for(j=0;j<3;j++) out<<dev[3*i+j]<<" ";
    for(j=0;j<3;j++) out<<dev_sq[3*i+j]<<" ";
    out<<std::endl;
  }
  out.close();
};


void GeneralSimulator::expansion(double T,double *scale) {
  double coeff;

  scale[0] = 1.0;
  coeff = std::stod(params->parameters["LinearThermalExpansionX"]);
  scale[0] += coeff*T;
  coeff = std::stod(params->parameters["QuadraticThermalExpansionX"]);
  scale[0] += coeff*T*T;

  scale[1] = 1.0;
  coeff = std::stod(params->parameters["LinearThermalExpansionY"]);
  scale[1] += coeff*T;
  coeff = std::stod(params->parameters["QuadraticThermalExpansionY"]);
  scale[1] += coeff*T*T;

  scale[2] = 1.0;
  coeff = std::stod(params->parameters["LinearThermalExpansionZ"]);
  scale[2] += coeff*T;
  coeff = std::stod(params->parameters["QuadraticThermalExpansionZ"]);
  scale[2] += coeff*T*T;

  //std::cout<<scale[0]<<" "<<scale[1]<<" "<<scale[2]<<std::endl;
};

void GeneralSimulator::make_path(std::vector<std::string> knot_list) {
  int i,j;
  int nknots = knot_list.size();
  // no way around it- have to store all the knots


  //std::vector<double> x(3*natoms,0.);

  std::vector<double> xs(nknots,0.), ys(nknots,0.), zs(nknots,0.);
  std::vector<double> r(nknots,0.);

  double *knots = new double[(const int)(3*nlocal*nknots)];
  double *lr = new double[(const int)(nknots)];

  spline xspl,yspl,zspl;
  double dx;

  // run through knots, and make spline
  load_config(knot_list[0],glo_data);
  for (i=0;i<3*nlocal;i++) loc_data[i] = glo_data[3*loc_id[i/3]+i%3];
  for (i=0;i<3*nlocal;i++) knots[i] = loc_data[i];

  for (int knot=1; knot<nknots; knot++) {
    load_config(knot_list[knot],glo_data);
    for (i=0;i<3*nlocal;i++) loc_data[i] = glo_data[3*loc_id[i/3]+i%3];
    for(i=0;i<3*nlocal;i++) loc_data[i]-=knots[i];
    pbc.wrap(loc_data);
    for(i=0;i<3*nlocal;i++) knots[i+knot*3*nlocal] = loc_data[i]+knots[i];
  }

  for(int knot=0;knot<nknots;knot++) {
    lr[knot] = 0.;
    for(i=0;i<3*nlocal;i++) {
      dx = knots[i+knot*3*nlocal]-knots[i];
      lr[knot] += dx*dx;
    }
  }
  MPI_Barrier(*comm);

  MPI_Allreduce(lr,&r[0],nknots,MPI_DOUBLE,MPI_SUM,*comm);

  for(int knot=0;knot<nknots;knot++) {
    r[knot] = sqrt(r[knot]/r[nknots-1]);
    //if(tag==0) std::cout<<local_rank<<" "<<knot<<" "<<r[knot]<<std::endl;
  }

  for(i=0; i<nlocal; i++) {
    for(int knot=0;knot<nknots;knot++) {
      xs[knot] = knots[3*nlocal*knot + 3*i + 0];
      ys[knot] = knots[3*nlocal*knot + 3*i + 1];
      zs[knot] = knots[3*nlocal*knot + 3*i + 2];
    }

    xspl.set_points(r,xs);
    pathway.push_back(xspl);

    yspl.set_points(r,ys);
    pathway.push_back(yspl);

    zspl.set_points(r,zs);
    pathway.push_back(zspl);

  }
  delete [] knots; // clear memory
  delete [] lr;

};

void GeneralSimulator::fill_path(double r,int der, std::vector<double> &vec){
  for(int i=0; i<3*natoms; i++) glo_data[i] = 0.0;
  for(int i=0; i<3*natoms; i++) vec[i]=0.0;

  for(int i=0; i<3*nlocal; i++) glo_data[3*loc_id[i/3]+i%3] = pathway[i].deriv(der,r);
  MPI_Barrier(*comm);

  MPI_Allreduce(&glo_data[0],&vec[0],3*natoms,MPI_DOUBLE,MPI_SUM,*comm);
};
