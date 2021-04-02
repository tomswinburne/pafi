#ifndef CGA_H
#define CGA_H

#include "GenericGatherer.hpp"

class TI_Gatherer : public GenericGatherer {
public:
  TI_Gatherer(Parser &p, std::vector<double> pathway_r, int _nW, int di, int _rank) :
    GenericGatherer(p, pathway_r, _nW, di, _rank){};

  // generated whenever a temperature cycle is complete
  void screen_output_header(bool end=true) {
    GenericGatherer::screen_output_header(false); // false : no endl;
    if(rank>0) return; // for flexibility
    std::cout<<std::setw(fw)<<" Opt"<<std::endl;
  };

  // generated whenever a reaction coordinate cycle is complete
  void screen_output_line(bool end=true) {
    GenericGatherer::screen_output_line(false); // false : no endl;
    if(rank>0) return; // for flexibility
    std::cout<<std::setw(fw)<<" 0.0"<<std::endl;
  };

};

#endif
