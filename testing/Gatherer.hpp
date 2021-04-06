#ifndef CGA_H
#define CGA_H

#include "GenericGatherer.hpp"

class TI_Gatherer : public GenericGatherer {
public:
  TI_Gatherer(Parser &p, int _nW, int di, int _rank) :
    GenericGatherer(p, _nW, di, _rank){
      field_width=15;
    };

  // generated whenever a temperature cycle is complete
  void screen_output_header(bool end=true) {
    GenericGatherer::screen_output_header(false); // false : no endl;
    if(rank>0) return; // for flexibility
    std::cout<<std::setw(field_width)<<"Lambda"<<std::endl;
  };

  // generated whenever a reaction coordinate cycle is complete
  void screen_output_line(bool end=true) {
    GenericGatherer::screen_output_line(false); // false : no endl;
    if(rank>0) return; // for flexibility
    std::cout<<std::setw(field_width)<<params["Lambda"]<<std::endl;
  };

};

#endif
