#ifndef PBC_H
#define PBC_H

#include <Eigen/Dense> // For Supercell and Hessian

class MinImage {
public:
  std::array< std::array<double,3>, 3 > cell;
  std::array< std::array<double,3>, 3 > invcell;
  std::array<double,3> periodic;
  MinImage(){};

  void load(void *ptr) { // Overload to allow lammps initialization
    Eigen::MatrixXd c(3,3), invc(3,3); // temp as array has faster access

    LAMMPS_NS::LAMMPS *lmp = (LAMMPS_NS::LAMMPS *) ptr;

    double *boxlo = (double *) lammps_extract_global(lmp,(char *) "boxlo");
  	double *boxhi = (double *) lammps_extract_global(lmp,(char *) "boxhi");
    int *pv = (int *) lammps_extract_global(lmp,(char *) "periodicity");

    for(int i=0; i<3; i++) {
      periodic[i] = (double)pv[i];
      for(int j=0;j<3;j++) cell[i][j] = 0.;
      cell[i][i] = boxhi[i]-boxlo[i];
    }

    cell[0][1] = *((double *) lammps_extract_global(lmp,(char *) "xy"));
    cell[0][2] = *((double *) lammps_extract_global(lmp,(char *) "xz"));
    cell[1][2] = *((double *) lammps_extract_global(lmp,(char *) "yz"));

    for(int i=0; i<3; i++) for(int j=0; j<3; j++) c(i,j)=cell[i][j];
    invc=c.inverse();
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) invcell[i][j]=invc(i,j);
  };

  void minimumImageVector(double *dr) {
  	std::array<double, 3> ds;
  	//project to cell-space
  	for(int i=0; i<3; i++) {
  		ds[i]=0;
  		for(int j=0; j<3; j++) ds[i] += invcell[i][j]*dr[j];
  	}
  	//wrap if necessary (where allowed)
  	for(int i=0; i<3; i++)	ds[i] -= round(ds[i])*periodic[i];
  	//project to R-space
  	for(int i=0; i<3; i++) {
  		dr[i]=0;
  		for(int j=0; j<3; j++) dr[i] += cell[i][j]*ds[j];
  	}
  };

  void wrap(std::vector<double> &dx) {
    double t[3];
    for(int i=0;i<dx.size()/3;i++) {
      for(int j=0;j<3;j++) t[j] = dx[3*i+j];
      minimumImageVector(t);
      for(int j=0;j<3;j++) dx[3*i+j] = t[j];
    }
  };
};

#endif
