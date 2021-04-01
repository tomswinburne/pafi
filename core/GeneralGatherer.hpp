#ifndef UTILS_H
#define UTILS_H


#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <array>

#include "Parser.hpp"


class GenericGatherer {
  public:
    virtual int initialize(Parser &p, std::vector<double> pathway_r, int _nW, int di, int _rank) {
      rank = _rank;
      dump_index = di;
      parser = &p;
      dsize=0;
      ens_data==NULL;
      nWorkers = _nW;
      fw=18;

      sweep_order.clear();

      // determine the arrays:
      for(auto aa: parser->parameters) {
        sample_axes[aa.first] =  *(new std::vector<double>);
        if(std::get<2>(aa.second)>1) {
          double delta = std::get<1>(aa.second)-std::get<0>(aa.second);
          delta /= 1.0*(std::get<2>(aa.second)-1);
          for(int i=0;i<std::get<2>(aa.second);i++) {
            sample_axes[aa.first].push_back(std::get<0>(aa.second) + i*delta);
          }
        } else sample_axes[aa.first].push_back(std::get<0>(aa.second));
        sweep_order.push_back(std::make_pair(aa.first,0));
      }

      if(sample_axes.find("Temperature")==sample_axes.end()) {
        std::cout<<"Must declare a Temperature parameter!!!"<<std::endl;
        return 0;
      }
      if(sample_axes.find("ReactionCoordinate")==sample_axes.end()) {
        std::cout<<"Must declare a ReactionCoordinate parameter!!!"<<std::endl;
        return 0;
      }

      // special option
      if((!parser->spline_path) or parser->match_planes) {
        sample_axes["ReactionCoordinate"].clear();
        for(auto r: pathway_r) if(r>=0.0 && r<=1.0) {
          sample_axes["ReactionCoordinate"].push_back(r);
        }
      }

      raw_dump_file = parser->dump_dir+"/raw_data_output_"+std::to_string(dump_index);
      if(rank==0) {
        if(raw.is_open()) raw.close();
        raw.open(raw_dump_file.c_str(),std::ofstream::out);
        if(raw.is_open()) return 1;
        return 0;
      } else return 1;
    };

    virtual std::map<std::string,double> params() {
      std::map<std::string,double> h;
      dev_file = parser->dump_dir+"/deviation_";
      for(auto s: sweep_order) {
        h[s.first] = sample_axes[s.first][s.second];
        dev_file += s.first+"_"+std::to_string(sample_axes[s.first][s.second]);
      }
      dev_file +="_"+std::to_string(dump_index);
      return h;
    }

    virtual void prepare(std::map<std::string,double> &sim_results) {
      if(!raw.is_open() or rank>0) return;
      int i=0,j=0;
      raw<<"# ";
      for(auto par: params()) raw<<i++<<": "<<par.first<<" ";
      for(auto res: sim_results) {
        raw<<i++<<": "<<res.first<<"  ";
        ens_results[res.first] = std::make_pair(0.,0.);
      }
      raw<<std::endl;

      dsize = sim_results.size();
      if(rank==0) {
        if(ens_data==NULL) ens_data = new double[dsize*2+1];
        // ensemble average
        for(j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
      }
    };

    virtual int collate(int *valid, double *all_data) {

      if(!raw.is_open() or rank>0) return 0;

      int total_valid=0,i,j;

      // raw output
      for(auto par: params()) raw<<par.second<<" ";
      for(i=0;i<nWorkers*dsize;i++) raw<<all_data[i]<<" ";
      raw<<std::endl;

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


      // put in results

      i=0;
      for(auto &res : ens_results) {
        res.second = std::make_pair(ens_data[i],ens_data[dsize+i]);
        i++;
      }


      /*
      data_log.push_back(r);
      if(end) for(auto f: log_fields) if(results.find(f.first)!=results.end()) {
        data_log.push_back(results[f.first]);
      } else if(local_rank==0) {
        std::cout<<"Field "<<f.first<<" not found!"<<std::endl;
      }
      */


      return total_valid;
    };

    virtual void next() {

      screen_output_line();

      for(auto &s : sweep_order) {
        s.second++;
        if(s.second % sample_axes[s.first].size() == 0) {
          // axis complete - trigger iteration...
          s.second = 0;
          if(rank==0) {
            std::cout<<"\n"<<s.first<<" : ";
            std::cout<<*(sample_axes[s.first].begin())<<" -> ";
            std::cout<<*std::prev(sample_axes[s.first].end());
            std::cout<<" iteration complete"<<std::endl;
          }

        } else {
          if(s.first=="Temperature" ) screen_output_header();
          break;
        }

      }
      // wipe ens_data if master node
      if(rank==0) for(int j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
    };

    virtual void screen_output_header() {
      if(rank>0) return;
      std::cout<<"\nStarting T="<<params()["Temperature"]<<"K run\n";
      std::cout<<std::setw(35)<<"r";
      std::cout<<std::setw(fw)<<"av(<Tpre>)";
      std::cout<<std::setw(fw)<<"av(<Tpost>)";
      std::cout<<std::setw(fw)<<"av(<dF/dr>)";
      std::cout<<std::setw(fw)<<"err(<dF/dr>)";
      std::cout<<std::setw(fw)<<"av(|(<X>-U).N|)";
      std::cout<<std::setw(fw)<<"av(<N_true>.N)";
      std::cout<<std::setw(fw)<<"Max Jump";
      std::cout<<std::setw(fw)<<"% Valid";
      std::cout<<std::endl;
    };

    virtual void screen_output_line() {
      if(rank>0) return;
      // screen output
      std::cout<<std::setw(35)<<params()["ReactionCoordinate"];//"r"
      std::cout<<std::setw(fw)<<ens_results["preT"].first;//"av(<Tpre>)"
      std::cout<<std::setw(fw)<<ens_results["postT"].first;//"av(<Tpost>)"
      std::cout<<std::setw(fw)<<ens_results["aveF"].first;//"av(<dF/dr>)"
      std::cout<<std::setw(fw)<<ens_results["aveF"].second;//"err(<dF/dr>)"
      std::cout<<std::setw(fw)<<ens_results["dXTangent"].first;//"av(|<X>-U).(dU/dr)|)"
      std::cout<<std::setw(fw)<<ens_results["avePsi"].first;//"av(Psi)"
      std::cout<<std::setw(fw)<<ens_results["MaxJump"].first;// max jump
      std::cout<<std::setw(fw)<<ens_results["Valid"].first*100.0;// ratio of jumps
      std::cout<<std::endl;
    };



    bool finished() {
      int isum=0;
      for(auto &s : sweep_order) isum += s.second;
      if(isum==0) return true;
      return false;
    };

    void close() {
      if(raw.is_open()) raw.close();
    };

  int dsize,nWorkers,dump_index,rank,fw;
  Parser *parser;
  std::ofstream raw;
  double *ens_data;
  std::map<std::string,std::vector<double>> sample_axes;
  std::map<std::string,std::pair<double,double>> ens_results;
  std::vector<std::pair<std::string,int>> sweep_order;
  std::string raw_dump_file, dump_suffix, dev_file;
};
#endif
