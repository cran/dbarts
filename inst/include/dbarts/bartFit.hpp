#ifndef DBARTS_BART_FIT_HPP
#define DBARTS_BART_FIT_HPP

#include <cstddef> // size_t
#include "cstdint.hpp" // uint32_t

#include <external/thread.h>

#include "types.hpp"

#include "control.hpp"
#include "data.hpp"
#include "model.hpp"
#include "scratch.hpp"
#include "state.hpp"

extern "C" {
  struct _ext_htm_manager_t;
  typedef struct _ext_htm_manager_t* ext_htm_manager_t;
}

namespace dbarts {
  struct Results;
  struct SharedScratch;
  
  struct BARTFit {
    Control control; // top three are passed in from elsewhere
    Model model;
    Data data;
    
    SharedScratch sharedScratch;
    ChainScratch* chainScratch;
    State* state;
    
    double runningTime;
    std::size_t currentNumSamples;
    std::size_t currentSampleNum;
    
    ext_htm_manager_t threadManager;
    
    BARTFit(Control control, Model model, Data data);
    ~BARTFit();
    
    void setRNGState(const void* const* uniformState, const void* const* normalState);
    
    Results* runSampler();
    Results* runSampler(std::size_t numBurnIn, std::size_t numSamples);
    void runSampler(std::size_t numBurnIn, Results* results);
    
    
    void predict(const double* x_test, std::size_t numTestObservations, const double* testOffset, double* result) const;
    // settors simply replace local pointers to variables. dimensions much match
    // update modifies the local copy (which may belong to someone else)
    void setResponse(const double* newResponse); 
    void setOffset(const double* newOffset);
    
    // predictor changes will return false if the new covariates would leave the sampler in an invalid state
    // (i.e. with an empty terminal node); the update functions auto-revert to the previous while set does not
    bool setPredictor(const double* newPredictor);
    bool updatePredictor(const double* newPredictor, std::size_t column); // false if same, but reverts on own
    bool updatePredictors(const double* newPredictor, const std::size_t* columns, std::size_t numColumns); 
    
    void setTestPredictor(const double* newTestPredictor, std::size_t numTestObservations);
    void setTestOffset(const double* newTestOffset);
    void setTestPredictorAndOffset(const double* newTestPredictor, const double* newTestOffset, std::size_t numTestObservations);
    
    void updateTestPredictor(const double* newTestPredictor, std::size_t column);
    void updateTestPredictors(const double* newTestPredictor, const std::size_t* columns, std::size_t numColumns);
    
    void sampleTreesFromPrior();
    
    void rebuildScratchFromState();
    
    void printTrees(const std::size_t* chains, std::size_t numChains,
                    const std::size_t* samples, std::size_t numSamples,
                    const std::size_t* indices, std::size_t numIndices) const;
    
    // this assumes that the new data has as many predictors as the old, and that they correspond to each other;
    // it'll attempt to map cut points from the old to the new, and prune any trees that may have been left in an
    // invalid state
    void setData(const Data& data);
    // the new control must have the same number of chains as the previous or else prob seg fault
    void setControl(const Control& control);
    void setModel(const Model& model);
    
    bool saveToFile(const char* fileName) const;
    static BARTFit* loadFromFile(const char* fileName);
  };
} // namespace dbarts

#endif // DBARTS_BART_FIT_HPP
