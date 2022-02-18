#include "AdaptiveFilter.h"

namespace gplib
  {

    AdaptiveFilter::AdaptiveFilter(const int inputsize, const int outputsize) :
      FilterOutput(outputsize), Epsilon(outputsize), inputlength(inputsize),
          outputlength(outputsize)
      {
      }
    AdaptiveFilter::~AdaptiveFilter()
      {
      }
  }
