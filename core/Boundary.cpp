#include "Boundary.hpp"

void MinImage::load(std::array<double,9> lcell) {
  determinant = 0.0;
  for(int i=0; i<9; i++) cell[i/3][i%3] = 0.;

  for(int i=0; i<3; i++) {
    cell[i][i] = lcell[i];
    cell[int(i==2)][1+int(i>0)] = lcell[3+i];
    periodic[i] = lcell[6+i];
  }
  
  determinant  = cell[0][0] * (cell[1][1]*cell[2][2] - cell[1][2] * cell[2][1]);
  determinant -= cell[0][1] * (cell[1][0]*cell[2][2] - cell[1][2] * cell[2][0]);
  determinant += cell[0][2] * (cell[1][0]*cell[2][1] - cell[1][1] * cell[2][0]);
  
  invcell[0][0] = (1.0/determinant) * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1]);
  invcell[1][0] = (1.0/determinant) * (cell[1][2] * cell[2][0] - cell[1][0] * cell[2][2]);
  invcell[2][0] = (1.0/determinant) * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0]);
  invcell[0][1] = (1.0/determinant) * (cell[0][2] * cell[2][1] - cell[0][1] * cell[2][2]);
  invcell[1][1] = (1.0/determinant) * (cell[0][0] * cell[2][2] - cell[0][2] * cell[2][0]);
  invcell[2][1] = (1.0/determinant) * (cell[0][1] * cell[2][0] - cell[0][0] * cell[2][1]);
  invcell[0][2] = (1.0/determinant) * (cell[0][1] * cell[1][2] - cell[0][2] * cell[1][1]);
  invcell[1][2] = (1.0/determinant) * (cell[0][2] * cell[1][0] - cell[0][0] * cell[1][2]);
  invcell[2][2] = (1.0/determinant) * (cell[0][0] * cell[1][1] - cell[0][1] * cell[1][0]);
};

/*
void MinImage::load(std::array<double,9> lcell) {
  Eigen::MatrixXd c(3,3), invc(3,3); // temp as array has faster access
  for(int i=0; i<9; i++) {
    cell[i/3][i%3] = 0.;
    invcell[i/3][i%3] = 0.;
    c(i/3,i%3) = 0.0;
    invc(i/3,i%3) = 0.0;
  }
  for(int i=0; i<3; i++) {
    cell[i][i] = lcell[i];
    cell[int(i==2)][1+int(i>0)] = lcell[3+i];
    periodic[i] = lcell[6+i];
  }
  for(int i=0; i<9; i++) c(i/3,i%3)=cell[i/3][i%3];
  invc=c.inverse();
  for(int i=0; i<9; i++) invcell[i/3][i%3]=invc(i/3,i%3);
};
*/

void MinImage::minimumImageVector(double *dr) {
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

void MinImage::wrap(std::vector<double> &dx) {
  double t[3];
  for(int i=0;i<dx.size()/3;i++) {
    for(int j=0;j<3;j++) t[j] = dx[3*i+j];
    minimumImageVector(t);
    for(int j=0;j<3;j++) dx[3*i+j] = t[j];
  }
};
