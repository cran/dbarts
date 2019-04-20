// stuff that can be figured out from Control, Model, and Data
#ifndef DBARTS_SCRATCH_HPP
#define DBARTS_SCRATCH_HPP

#include <cstddef> // size_t
#include <dbarts/types.hpp>

namespace dbarts {
  struct ScaleFactor { double min, max, range; };
  
  struct SharedScratch {
    const double* yRescaled;
    const xint_t* x; // x mapped to integers
    const xint_t* xt_test; // xt for the test obs
    
    ScaleFactor dataScale;
  };
  struct ChainScratch {
    double* treeY;
    double* probitLatents;
    
    double* totalFits;     // numObs
    double* totalTestFits; // numTestObs
    
    std::size_t taskId;
  };
} // namespace dbarts

#endif // DBARTS_SCRATCH_HPP

