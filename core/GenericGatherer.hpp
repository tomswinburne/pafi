#ifndef UTILS_H
#define UTILS_H


#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <array>

#include "ConstantsTypes.hpp"
#include "Parser.hpp"


class GenericGatherer {
public:
GenericGatherer(Parser &p, int _nW, int di, int _rank) {
  initialized = 1;
  rank = _rank;
  dump_index = di;
  parser = &p;
  dsize=0;
  ens_data=NULL;
  data=NULL;
  all_data=NULL;
  nWorkers = _nW;
  field_width=18;

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
    initialized = 0;
  }
  if(sample_axes.find("ReactionCoordinate")==sample_axes.end()) {
    std::cout<<"Must declare a ReactionCoordinate parameter!!!"<<std::endl;
    initialized = 0;
  }

  raw_dump_file = parser->dump_dir+"/raw_data_output_"+std::to_string(dump_index);

  if(rank==0 and initialized==1) {
    if(raw.is_open()) raw.close();
    raw.open(raw_dump_file.c_str(),std::ofstream::out);
    if(!raw.is_open()) {
      std::cout<<"Could not open raw dump file! Does the directory";
      std::cout<<parser->dump_dir<<" Exist? EXIT"<<std::endl;
      initialized = 0;
    }
  }
  set_params();
};

void special_r_overwrite(std::vector<double> pathway_r) {
  // special option
  if((!parser->spline_path) or parser->match_planes) {
    sample_axes["ReactionCoordinate"].clear();
    for(auto r: pathway_r) if(r>=0.0 && r<=1.0) {
      sample_axes["ReactionCoordinate"].push_back(r);
    }
  }
};

virtual void set_params() {
  dev_file = parser->dump_dir+"/deviation_";
  for(auto s: sweep_order) {
    params[s.first] = sample_axes[s.first][s.second];
    dev_file += s.first+"_"+std::to_string(sample_axes[s.first][s.second]);
  }
  dev_file +="_"+std::to_string(dump_index);
}

virtual void prepare(Holder &sim_results) {

  int j=0, i=0;

  if(dsize==0) {

    dsize = sim_results.size(); // constant for a given simulation

    data = new double[dsize];

    if (rank==0 and raw.is_open()) {
      // data holders

      ens_data = new double[dsize*2+1];
      for(j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;

      all_data = new double[dsize*nWorkers];
      for(j=0;j<dsize*nWorkers;j++) all_data[j] = 0.0;

      for(auto res: sim_results)
        ens_results[res.first] = *(new std::pair<double,double>);

      raw<<"# ";
      i = -params.size();
      for(auto par: params) raw<<i++<<": "<<par.first<<" ";
      for(auto res: sim_results) raw<<i++<<": "<<res.first<<"  ";
      raw<<std::endl;
    }

  }

  i=0; for(auto res: sim_results) data[i++] = res.second; // fill

};

virtual int collate(int *valid) {

  if(!raw.is_open() or rank>0) return 0;

  int total_valid=0,i,j;

  // raw output
  for(auto par: params) raw<<par.second<<" ";
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
  i=0; for(auto res = ens_results.begin(); res!=ens_results.end(); res++, i++)
    res->second = std::make_pair(ens_data[i],ens_data[dsize+i]);
  all_ens_results.insert(std::make_pair(params,ens_results));

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
        std::cout<<s.first<<" : ";
        std::cout<<*(sample_axes[s.first].begin())<<" -> ";
        std::cout<<*std::prev(sample_axes[s.first].end());
        std::cout<<" complete"<<std::endl;
      }

    } else {
      if(s.first=="Temperature" ) screen_output_header();
      break;
    }

  }
  set_params();
  // wipe ens_data if master node
  if(rank==0) {
    for(int j=0;j<2*dsize+1;j++) ens_data[j] = 0.0;
  }
};

virtual void screen_output_header(bool end=true) {
  if(rank>0) return;
  std::cout<<"\nStarting T="<<params["Temperature"]<<"K run\n";
  std::cout<<std::setw(35)<<"r";
  std::cout<<std::setw(field_width)<<"av(<Tpre>)";
  std::cout<<std::setw(field_width)<<"av(<Tpost>)";
  std::cout<<std::setw(field_width)<<"av(<dF/dr>)";
  std::cout<<std::setw(field_width)<<"err(<dF/dr>)";
  std::cout<<std::setw(field_width)<<"av(|(<X>-U).N|)";
  std::cout<<std::setw(field_width)<<"av(<N_true>.N)";
  std::cout<<std::setw(field_width)<<"Max Jump";
  std::cout<<std::setw(field_width)<<"% Valid";
  if(end) std::cout<<std::endl;
};

virtual void screen_output_line(bool end=true) {
  if(rank>0) return;
  // screen output
  std::cout<<std::setw(35)<<params["ReactionCoordinate"];//"r"
  std::cout<<std::setw(field_width)<<ens_results["preT"].first;//"av(<Tpre>)"
  std::cout<<std::setw(field_width)<<ens_results["postT"].first;//"av(<Tpost>)"
  std::cout<<std::setw(field_width)<<ens_results["aveF"].first;//"av(<dF/dr>)"
  std::cout<<std::setw(field_width)<<ens_results["aveF"].second;//"err(<dF/dr>)"
  std::cout<<std::setw(field_width)<<ens_results["dXTangent"].first;//"av(|<X>-U).(dU/dr)|)"
  std::cout<<std::setw(field_width)<<ens_results["avePsi"].first;//"av(Psi)"
  std::cout<<std::setw(field_width)<<ens_results["MaxJump"].first;// max jump
  std::cout<<std::setw(field_width)<<ens_results["Valid"].first*100.0;// ratio of jumps
  if(end) std::cout<<std::endl;
};

bool finished() {
  int isum=0;
  for(auto &s : sweep_order) isum += s.second;
  if(isum==0) return true;
  return false;
};

void close() {
  if(raw.is_open()) raw.close();
  // perform any integrations here....
};

int dsize,nWorkers,dump_index,rank,field_width,initialized;
Holder params;
Parser *parser;
std::ofstream raw;
double *ens_data, *data, *all_data;
std::map<std::string,std::vector<double>> sample_axes;
std::vector<std::pair<std::string,int>> sweep_order;
PairHolder<double,double> ens_results;
std::map<Holder,PairHolder<double,double>> all_ens_results;
std::string raw_dump_file, dump_suffix, dev_file;
};
#endif
