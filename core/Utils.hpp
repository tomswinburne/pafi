#ifndef UTILS_H
#define UTILS_H


#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <array>

//#include <Eigen/Dense> // For Supercell and Hessian

class DataGatherer {
  public:
    DataGatherer(){
      dsize=0;
    }

    int initialize(Parser &p, std::string dump_file,int _nWorkers) {
      if(raw.is_open()) raw.close();
      nWorkers = _nWorkers;
      params = &p;
      raw.open(dump_file.c_str(),std::ofstream::out);
      if(raw.is_open()) return 1;
      return 0;
    };

    void prepare(std::map<std::string,double> &results) {
      if(!raw.is_open()) return;
      int i=1,j=0;
      std::vector<double> blank;

      raw<<"# 0: r ";
      for(auto res: results) {
        raw<<i++<<": "<<res.first<<"  ";
        all_results[res.first] = blank;
      }
      all_results["sample_r"] = blank;
      raw<<std::endl;

      dsize = results.size();
      all_data = new double[dsize*nWorkers];
      ens_data = new double[dsize*2+1];


      // ensemble average
      for(j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
    };

    int ensemble(double r, int *valid) {
      if(!raw.is_open()) return 0;

      int total_valid=0,i,j;

      // raw output
      raw<<r<<" ";
      for(i=0;i<nWorkers*dsize;i++) raw<<all_data[i]<<" ";
      raw<<std::endl;

      // all results
      j=0;
      all_results["sample_r"].push_back(r);
      for(auto &res: all_results) {
        for(i=0;i<nWorkers;i++) res.second.push_back(all_data[i*dsize+j]);
        j++;
      }

      // unnormalize for sum
      for(j=0;j<dsize;j++) {
        ens_data[j+dsize] *= ens_data[2*dsize];
        ens_data[j+dsize] += ens_data[j] * ens_data[j];
      }
      for(j=0;j<2*dsize;j++) ens_data[j] *= ens_data[2*dsize];

      // add more data
      for(i=0;i<nWorkers;i++) if(valid[i]==1) {
        total_valid++;
        ens_data[2*dsize] += 1.0;
        for(j=0;j<dsize;j++) {
          ens_data[j]+=all_data[i*dsize+j];
          ens_data[j+dsize] += all_data[i*dsize+j] * all_data[i*dsize+j];
        }
      }

      // renormalize for access
      if(ens_data[2*dsize]>0.5) {
        for(j=0;j<2*dsize;j++) ens_data[j] /= ens_data[2*dsize];
        for(j=0;j<dsize;j++) {
          ens_data[j+dsize] -= ens_data[j] * ens_data[j];
          ens_data[j+dsize] /= ens_data[2*dsize]; // N^2 for aves-of-aves
        }
      }
      return total_valid;
    };

    void next() {
      // wipe ens_data
      for(int j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
    };

    void close() {
      if(raw.is_open()) raw.close();
    };


  Parser *params;
  std::ofstream raw;
  double *all_data,*ens_data;
  int dsize,nWorkers;
  std::map<std::string,std::vector<double>> all_results;
};
#endif
