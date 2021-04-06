#ifndef CGA_H
#define CGA_H

#include "GenericGatherer.hpp"

class StandardGatherer : public GenericGatherer {
public:
  StandardGatherer(Parser &p, int _nW, int di, int _rank) :
    GenericGatherer(p, _nW, di, _rank){};

  // generated whenever a temperature cycle is complete
  void screen_output_header(bool end=true) {
    /*
      Header line for screen output
      GenericGatherer::screen_output_header only outputs for rank==0,
      but this function will output for any rank
      boolean argument : linebreak at end of output.
      If false, other arguments can be appended. See TI example
    */
    GenericGatherer::screen_output_header(true); // false : no endl;
  };

  // generated whenever a reaction coordinate cycle is complete
  void screen_output_line(bool end=true) {
    /*
      Header line for screen output of a single sampling run
      GenericGatherer::screen_output_line only outputs for rank==0,
      but this function will output for any rank
      boolean argument : linebreak at end of output.
      If false, other arguments can be appended. See TI example
     */
    GenericGatherer::screen_output_line(true);
  };

};

#endif
