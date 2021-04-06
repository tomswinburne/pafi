#ifndef CGA_H
#define CGA_H

#include "GenericGatherer.hpp"

class LambdaGatherer : public GenericGatherer {
public:
  LambdaGatherer(Parser &p, int _nW, int di, int _rank) :
    GenericGatherer(p, _nW, di, _rank){
      field_width=15; // narrower print out (more fields)
    };

  void screen_output_header(bool end=true) {
    /*
      Header line for screen output
      GenericGatherer::screen_output_header only outputs for rank==0,
      but this function will output for any rank
      boolean argument : linebreak at end of output.
      If false, other arguments can be appended.
    */
    GenericGatherer::screen_output_header(false);
    if(rank>0) return; // only print for rank==0
    std::cout<<std::setw(field_width)<<"Lambda"<<std::endl;
  };

  void screen_output_line(bool end=true) {
    /*
      Header line for screen output of a single sampling run
      GenericGatherer::screen_output_line only outputs for rank==0,
      but this function will output for any rank
      boolean argument : linebreak at end of output.
      If false, other arguments can be appended.
     */
    GenericGatherer::screen_output_line(false);
    if(rank>0) return; // only print for rank==0
    std::cout<<std::setw(field_width)<<params["Lambda"]<<std::endl;
  };
};

#endif
