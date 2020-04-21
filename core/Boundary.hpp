#ifndef PBC_H
#define PBC_H


#include <math.h>
#include <vector>
#include <string>
#include <array>

//#include <Eigen/Dense> // For Supercell and Hessian



class MinImage {
public:
  double volume;
  std::array< std::array<double,3>, 3 > cell;
  std::array< std::array<double,3>, 3 > invcell;
  std::array<double,3> periodic;
  MinImage(){};
  void load(std::array<double,9> lcell);
  void minimumImageVector(double *dr);
  void wrap(double *dx,int dxsize);
};
#endif
