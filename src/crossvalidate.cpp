#include "config.hpp"
#include "crossvalidate.hpp"

#include <algorithm> // sort
#include <cstddef> // size_t
#include <cmath>   // sqrt, floor

#include <external/alloca.h>
#include <external/io.h>
#include <external/linearAlgebra.h>
#include <external/random.h>
#include <external/stats.h>
#include <external/thread.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>

#define INVALID_INDEX static_cast<size_t>(-1)

using std::size_t;

namespace {
  using namespace dbarts;
  using namespace xval;
  
  void permuteIndexArray(ext_rng* restrict generator, size_t* restrict indices, size_t length);
  
  void allocateDataStorage(const Data& origData, Data& repData, size_t numTrainingSamples, size_t numTestSamples);
  void allocateModelStorage(const Model& origModel, Model& repModel);
  void divideData(const Data& restrict origData, Data& restrict repData, double* restrict y_test,
                  size_t numTrainingObservations, size_t numTestObservations,
                  ext_rng* restrict generator, size_t* restrict permutation);

  void freeDataStorage(Data& repData);
  void freeModelStorage(Model& repModel);

  struct CellParameters {
    size_t numTrees;
    double k;
    double power;
    double base;
    uint8_t status;
  };
  
  void updateFitForCell(BARTFit& fit, Control& repControl, Model& repModel, const CellParameters& parameters, bool verbose);
  
  struct SharedData {
    const Control& control;
    const Model& model;
    const Data& data;
    
    size_t numInitialBurnIn;
    size_t numContextShiftBurnIn;
    size_t numRepBurnIn;
    
    size_t numTrainingObservations;
    size_t numTestObservations;
    
    const LossFunctorDefinition& lossFunctorDef;
    
    size_t numReps;
    const CellParameters* parameters;
    
    size_t* printedCells;
  };
  
  struct ThreadData {
    SharedData* shared;
    ext_rng* rng;
    
    size_t repCellOffset;
    size_t numRepCells;
    size_t fittingCell;
    
    double* results;
  };
}

extern "C" { static void crossvalidationTask(void* data); }
extern "C" { static void printInfo(void** data, size_t numThreads); }

namespace dbarts { namespace xval {
    void crossvalidate(const Control& origControl, const Model& origModel, const Data& origData,
                       double testSampleProp, size_t numReps,
                       size_t numInitialBurnIn, size_t numContextShiftBurnIn, size_t numRepBurnIn,
                       const LossFunctorDefinition& lossFunctorDef, size_t numThreads,
                       const std::size_t* nTrees, size_t numNTrees, const double* k, size_t numKs,
                       const double* power, size_t numPowers, const double* base, size_t numBases,
                       double* results)

  {
    size_t numTestObservations     = origData.numObservations -
                                     static_cast<size_t>(std::floor(static_cast<double>(origData.numObservations) * testSampleProp + 0.5));
    size_t numTrainingObservations = origData.numObservations - numTestObservations;
    
    if (origControl.verbose) {
      ext_printf("starting crossvalidation with %.2f%% test obs, %lu replications, %lu/%lu test/training obs\n",
                 100.0 * testSampleProp, numReps, numTestObservations, numTrainingObservations);
      ext_printf("  %lu tree par(s), %lu k par(s), %lu power par(s), %lu base par(s)\n",
                 numNTrees, numKs, numPowers, numBases);
      ext_printf("  results of type: %s\n", lossFunctorDef.displayString);
      ext_printf("  num samp: %lu, num reps: %lu\n", origControl.defaultNumSamples, numReps);
      ext_printf("  burn in: %lu first, %lu shift, %lu rep\n\n", numInitialBurnIn, numContextShiftBurnIn, numRepBurnIn);

      if (numThreads > 1) {
        ext_printf("  parameters for [thread num, cell number]; some cells may be split across multiple threads:\n\n");
      }
      ext_fflush_stdout();
    }
    
    size_t numRepCells = numNTrees * numKs * numPowers * numBases * numReps;
    if (numRepCells < numThreads) numThreads = numRepCells;
    
    size_t numCells = numRepCells / numReps;
    CellParameters* cellParameters = new CellParameters[numCells];
    
    size_t cellNumber = 0;
    for (size_t nIndex = 0; nIndex < numNTrees; ++nIndex) {
      for (size_t kIndex = 0; kIndex < numKs; ++kIndex) {
        for (size_t pIndex = 0; pIndex < numPowers; ++pIndex) {
          for (size_t bIndex = 0; bIndex < numBases; ++bIndex) {
            
            cellParameters[cellNumber].numTrees = nTrees[nIndex];
            cellParameters[cellNumber].k        = k[kIndex];
            cellParameters[cellNumber].power    = power[pIndex];
            cellParameters[cellNumber].base     = base[bIndex];
            // cellParameters[cellNumber].offset   = nIndex + numNTrees * (kIndex + numKs * (pIndex + bIndex * numPowers));
            cellNumber++;
          }
        }
      }
    }
    
    Control threadControl = origControl;
    threadControl.verbose    = origControl.verbose == true && numThreads == 1;
    threadControl.numThreads = numThreads;
    
    SharedData sharedData = { threadControl, origModel, origData,
                              numInitialBurnIn, numContextShiftBurnIn, numRepBurnIn,
                              numTrainingObservations, numTestObservations,
                              lossFunctorDef, numReps, cellParameters,
                              new size_t[threadControl.numThreads] };
    for (size_t i = 0; i < numThreads; ++i) sharedData.printedCells[i] = INVALID_INDEX;
    
    if (numThreads <= 1) {
      ext_rng* rng;
      
      // both BART sampler and xval algorithm can use native sampler
      if (threadControl.rng_algorithm == EXT_RNG_ALGORITHM_INVALID &&
          threadControl.rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID) {
        
        if ((rng = ext_rng_createDefault(true)) == NULL)
          ext_throwError("could not allocate rng");
      
      } else {
        if (ext_rng_createAndSeed(&rng, threadControl.rng_algorithm, threadControl.rng_standardNormal) != 0)
          ext_throwError("could not allocate rng");
      }
      
      ThreadData threadData = { &sharedData, rng, 0, numRepCells, INVALID_INDEX, results };
      
      crossvalidationTask(&threadData);
      
      ext_rng_destroy(rng);
    } else {
      if (threadControl.rng_algorithm == EXT_RNG_ALGORITHM_INVALID)
        threadControl.rng_algorithm = ext_rng_getDefaultAlgorithmType();
      if (threadControl.rng_standardNormal == EXT_RNG_STANDARD_NORMAL_INVALID)
        threadControl.rng_standardNormal = ext_rng_getDefaultStandardNormalType();
      
      ext_mt_manager_t threadManager;
      ext_mt_create(&threadManager, numThreads);
      
      size_t numRepCellsPerThread;
      size_t offByOneIndex;
      
      ext_mt_getNumThreadsForJob(threadManager, numRepCells, 0, NULL, &numRepCellsPerThread, &offByOneIndex);
      
      
      ThreadData* threadData = ext_stackAllocate(numThreads, ThreadData);
      void** threadDataPtrs  = ext_stackAllocate(numThreads, void*);
      for (size_t i = 0; i < offByOneIndex; ++i) {
        threadData[i].shared = &sharedData;
        if (ext_rng_createAndSeed(&threadData[i].rng, threadControl.rng_algorithm, threadControl.rng_standardNormal) != 0)
          ext_throwError("could not allocate rng");
        threadData[i].repCellOffset = i * numRepCellsPerThread;
        threadData[i].numRepCells   = numRepCellsPerThread;
        threadData[i].fittingCell   = INVALID_INDEX;
        threadData[i].results = results + i * numRepCellsPerThread * lossFunctorDef.numResults;
        threadDataPtrs[i] = threadData + i;
      }
      
      for (size_t i = offByOneIndex; i < numThreads; ++i) {
        threadData[i].shared = &sharedData;
        if (ext_rng_createAndSeed(&threadData[i].rng, threadControl.rng_algorithm, threadControl.rng_standardNormal) != 0)
          ext_throwError("could not allocate rng");
        threadData[i].repCellOffset = offByOneIndex * numRepCellsPerThread + (i - offByOneIndex) * (numRepCellsPerThread - 1);
        threadData[i].numRepCells   = numRepCellsPerThread - 1;
        threadData[i].fittingCell   = INVALID_INDEX;
        threadData[i].results = results + threadData[i].repCellOffset * lossFunctorDef.numResults;
        threadDataPtrs[i] = threadData + i;
      }
      
      if (origControl.verbose == true)
        ext_mt_runTasksWithInfo(threadManager, crossvalidationTask, threadDataPtrs, numThreads, 1, printInfo);
      else
        ext_mt_runTasks(threadManager, crossvalidationTask, threadDataPtrs, numThreads);
      
      for (size_t i = 0; i < numThreads; ++i)
        ext_rng_destroy(threadData[i].rng);
         
      ext_stackFree(threadDataPtrs);
      ext_stackFree(threadData);
        
      ext_mt_destroy(threadManager);
      
    }
    
    delete [] sharedData.printedCells;
    
    delete [] cellParameters;
    
  }
} }

extern "C" {
  using namespace dbarts;
  using namespace dbarts::xval;
  
  static void crossvalidationTask(void* v_data)
  {
    ThreadData& threadData(*static_cast<ThreadData*>(v_data));
    const SharedData& sharedData(*threadData.shared);
    
    const Control& origControl(sharedData.control);
    const Model& origModel(sharedData.model);
    const Data& origData(sharedData.data);
    
    size_t numTrainingObservations = sharedData.numTrainingObservations;
    size_t numTestObservations     = sharedData.numTestObservations;
    size_t numSamples              = origControl.defaultNumSamples;
    
    const LossFunctorDefinition& lfDef(sharedData.lossFunctorDef);
    
    LossFunctor* lf = lfDef.createFunctor(lfDef, numTestObservations, numSamples);
    
    double* suppliedY_test      = lfDef.y_testOffset      >= 0 ?
                                  *reinterpret_cast<double**>(reinterpret_cast<char*>(lf) + lfDef.y_testOffset) :
                                  NULL;
    double* suppliedTestSamples = sharedData.lossFunctorDef.testSamplesOffset >= 0 ?
                                  *reinterpret_cast<double**>(reinterpret_cast<char*>(lf) + lfDef.testSamplesOffset) :
                                  NULL;
    Results* samples =
      suppliedTestSamples == NULL ?
        new Results(numTrainingObservations, origData.numPredictors, numTestObservations, numSamples, 1) :
        new Results(numTrainingObservations, origData.numPredictors, numTestObservations, numSamples, 1,
                    new double[origControl.defaultNumSamples],
                    new double[numTrainingObservations * numSamples],
                    suppliedTestSamples,
                    new double[origData.numPredictors * numSamples]);
    
    Control repControl = origControl;
    repControl.numThreads = 1;
    bool verbose = repControl.verbose;
    repControl.verbose = false;
    
    BARTFit* fit = new BARTFit(repControl, origModel, origData);
    
    Data repData;
    Model repModel;
    
    allocateDataStorage(fit->data, repData, numTrainingObservations, numTestObservations);
    allocateModelStorage(fit->model, repModel);
    
    double* y_test = (suppliedY_test == NULL ? new double[numTestObservations] : suppliedY_test);
    
    size_t* permutation = new size_t[origData.numObservations];
    for (size_t i = 0; i < origData.numObservations; ++i) permutation[i] = i;
    
    size_t firstCell    = threadData.repCellOffset / sharedData.numReps;
    size_t firstCellRep = threadData.repCellOffset % sharedData.numReps;
    size_t lastCell     = (threadData.repCellOffset + threadData.numRepCells) / sharedData.numReps;
    size_t lastCellRep  = (threadData.repCellOffset + threadData.numRepCells) % sharedData.numReps;
    
    // ext_printf("running task with rep offset: %lu, num reps: %lu, num rep cells: %lu\n", threadData.repCellOffset, sharedData.numReps, threadData.numRepCells);
    // ext_printf("first cell: %lu, first cell rep: %lu, last cell: %lu, last cell rep: %lu\n", firstCell, firstCellRep, lastCell, lastCellRep);
    
    size_t resultIndex = 0;
    size_t numBurnIn = sharedData.numInitialBurnIn;
    
    // first and last cells are a bit of a mess, since there can be a lot of off-by-one stuff
    if (firstCellRep != 0) {
      // ext_printf("updating for cell %lu\n", firstCell);
      updateFitForCell(*fit, repControl, repModel, sharedData.parameters[firstCell], verbose);
      threadData.fittingCell = firstCell;
      
      for (size_t repIndex = firstCellRep; repIndex < sharedData.numReps; ++repIndex)
      {
        divideData(origData, repData, y_test, numTrainingObservations, numTestObservations, threadData.rng, permutation);
        fit->setData(repData);
        
        fit->runSampler(numBurnIn, samples);
                
        lfDef.calculateLoss(*lf, y_test, numTestObservations, samples->testSamples, numSamples, threadData.results + resultIndex);
        resultIndex += lfDef.numResults;
        
        numBurnIn = sharedData.numRepBurnIn;
      }
      
      ++firstCell;
      firstCellRep = 0;
      
      numBurnIn = sharedData.numContextShiftBurnIn;
    }
    
    for (size_t cellIndex = firstCell; cellIndex < lastCell; ++cellIndex) {
      // ext_printf("updating for cell %lu\n", cellIndex);
      updateFitForCell(*fit, repControl, repModel, sharedData.parameters[cellIndex], verbose);
      threadData.fittingCell = cellIndex;
      
      for (size_t repIndex = 0; repIndex < sharedData.numReps; ++repIndex)
      {
        divideData(origData, repData, y_test, numTrainingObservations, numTestObservations, threadData.rng, permutation);
        // ext_printf("  perm: %lu", permutation[0]);
        // for (size_t i = 1; i < 20; ++i) ext_printf(" %lu", permutation[i]);
        // ext_printf("\n");
        fit->setData(repData);
        
        fit->runSampler(numBurnIn, samples);
              
        lfDef.calculateLoss(*lf, y_test, numTestObservations, samples->testSamples, numSamples, threadData.results + resultIndex);
        resultIndex += lfDef.numResults;
        
        numBurnIn = sharedData.numRepBurnIn;
      }
      
      numBurnIn = sharedData.numContextShiftBurnIn;
    }
    
    if (lastCellRep != 0) {
      // ext_printf("updating for cell %lu\n", lastCell);
      updateFitForCell(*fit, repControl, repModel, sharedData.parameters[lastCell], verbose);
      threadData.fittingCell = lastCell;
      
      for (size_t repIndex = 0; repIndex < lastCellRep; ++repIndex)
      {
        divideData(origData, repData, y_test, numTrainingObservations, numTestObservations, threadData.rng, permutation);
        fit->setData(repData);
        
        fit->runSampler(numBurnIn, samples);
        
        lfDef.calculateLoss(*lf, y_test, numTestObservations, samples->testSamples, numSamples, threadData.results + resultIndex);
        resultIndex += lfDef.numResults;
        
        numBurnIn = sharedData.numRepBurnIn;
      }
    }
    threadData.fittingCell = lastCell + 1;
    
    
    delete [] permutation;
    
    if (suppliedY_test == NULL) delete [] y_test;
    
    delete fit;
    
    freeModelStorage(repModel);
    freeDataStorage(repData);
    
    if (suppliedTestSamples != NULL) samples->testSamples = NULL;
    delete samples;
    
    sharedData.lossFunctorDef.deleteFunctor(lf);
  }
}


namespace {
  using namespace dbarts;
  
  void permuteIndexArray(ext_rng* restrict generator, size_t* restrict indices, size_t length)
  {
    size_t temp, swapPos;
    for (size_t i = 0; i < length - 1; ++i) {
      swapPos = static_cast<size_t>(ext_rng_simulateUnsignedIntegerUniformInRange(generator, i, length));
      
      temp = indices[i];
      indices[i] = indices[swapPos];
      indices[swapPos] = temp;
    }
  }
  
  void divideData(const Data& restrict origData, Data& restrict repData, double* restrict y_test,
                  size_t numTrainingObservations, size_t numTestObservations,
                  ext_rng* restrict generator, size_t* restrict permutation)
  {
    size_t i, j, obsIndex;
    double* restrict y = const_cast<double* restrict>(repData.y);
    double* restrict x = const_cast<double* restrict>(repData.x);
    double* restrict x_test = const_cast<double* restrict>(repData.x_test);
    
    
    permuteIndexArray(generator, permutation, origData.numObservations);
    std::sort(permutation, permutation + numTestObservations);
    std::sort(permutation + numTestObservations, permutation + origData.numObservations);
    
    
    for (i = 0; i < numTestObservations; ++i) {
      obsIndex = *permutation++;
      y_test[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x_test[i + j * numTestObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
    }
    for (i = 0; i < numTrainingObservations; ++i) {
      obsIndex = *permutation++;
      y[i] = origData.y[obsIndex];
      for (j = 0; j < origData.numPredictors; ++j) {
        x[i + j * numTrainingObservations] = origData.x[obsIndex + j * origData.numObservations];
      }
    }
  }
  
  void allocateDataStorage(const Data& origData, Data& repData, size_t numTrainingObservations, size_t numTestObservations)
  {
    repData.y = new double[numTrainingObservations];
    repData.x = new double[numTrainingObservations * origData.numPredictors];
    repData.x_test = new double[numTestObservations * origData.numPredictors];
     
    repData.weights    = origData.weights != NULL ? new double[numTrainingObservations]  : NULL;
    repData.offset     = origData.offset  != NULL ? new double[numTrainingObservations]  : NULL;
    repData.testOffset = origData.offset  != NULL ? new double[numTestObservations] : NULL;
    
    repData.numObservations     = numTrainingObservations;
    repData.numPredictors       = origData.numPredictors;
    repData.numTestObservations = numTestObservations;
    repData.sigmaEstimate       = origData.sigmaEstimate;
    
    repData.variableTypes = origData.variableTypes;
    repData.maxNumCuts    = origData.maxNumCuts;
  }
  
  void allocateModelStorage(const Model& origModel, Model& repModel)
  {
    repModel.birthOrDeathProbability = origModel.birthOrDeathProbability;
    repModel.swapProbability = origModel.swapProbability;
    repModel.changeProbability = origModel.changeProbability;
    repModel.birthProbability = origModel.birthProbability;
    
    CGMPrior* repTreePrior = new CGMPrior();
    const CGMPrior* oldTreePrior = static_cast<CGMPrior*>(origModel.treePrior);
    repTreePrior->base = oldTreePrior->base;
    repTreePrior->power = oldTreePrior->power;
    
    repModel.treePrior = repTreePrior;
    
    
    NormalPrior* repNodePrior = new NormalPrior();
    const NormalPrior* oldNodePrior = static_cast<NormalPrior*>(origModel.muPrior);
    repNodePrior->precision = oldNodePrior->precision;
    
    repModel.muPrior = repNodePrior;
    
    
    ChiSquaredPrior* repResidPrior = new ChiSquaredPrior();
    const ChiSquaredPrior* oldResidPrior = static_cast<ChiSquaredPrior*>(origModel.sigmaSqPrior);
    repResidPrior->degreesOfFreedom = oldResidPrior->degreesOfFreedom;
    repResidPrior->scale = oldResidPrior->scale;
    
    repModel.sigmaSqPrior = repResidPrior;
  }
  
  void freeDataStorage(Data& repData)
  {
    delete [] repData.testOffset; repData.testOffset = NULL;
    delete [] repData.offset;     repData.offset     = NULL;
    delete [] repData.weights;    repData.weights    = NULL;
    
    delete [] repData.x_test; repData.x_test = NULL;
    delete [] repData.x;      repData.x      = NULL;
    delete [] repData.y;      repData.y      = NULL;
  }
  
  void freeModelStorage(Model& repModel)
  {
    delete repModel.sigmaSqPrior;
    delete repModel.muPrior;
    delete repModel.treePrior;
  }
}

extern "C" {
  void printInfo(void** vt_data, size_t numThreads)
  {
    SharedData& sharedData(*static_cast<ThreadData*>(vt_data[0])->shared);
    
    for (size_t i = 0; i < numThreads; ++i) {
      const ThreadData& threadData(*static_cast<const ThreadData*>(vt_data[i]));
      
      size_t fittingCell = threadData.fittingCell;
      size_t printedCell = sharedData.printedCells[i];
      
      
      if (fittingCell == INVALID_INDEX || (printedCell != INVALID_INDEX && printedCell >= fittingCell)) continue;
      
      size_t firstCell = threadData.repCellOffset / sharedData.numReps;
      size_t lastCell  = (threadData.repCellOffset + threadData.numRepCells) / sharedData.numReps;
      size_t lastCellRep  = (threadData.repCellOffset + threadData.numRepCells) % sharedData.numReps;
      if (lastCellRep == 0) --lastCell;
      
      size_t firstPrintCell = printedCell == INVALID_INDEX ? firstCell : printedCell + 1;
      size_t lastPrintCell = fittingCell < lastCell ? fittingCell : lastCell;
      
      if (firstPrintCell <= lastPrintCell) {
        for (size_t j = firstPrintCell; j <= lastPrintCell; ++j) {
          ext_printf("    [%lu, %lu] n.tree: %lu, k: %f, power: %f, base: %f\n",
                     i + 1, j + 1,
                     sharedData.parameters[j].numTrees, sharedData.parameters[j].k,
                     sharedData.parameters[j].power, sharedData.parameters[j].base);
        }
        sharedData.printedCells[i] = fittingCell;
      }
    }
  }
}

namespace {
  void updateFitForCell(BARTFit& fit, Control& repControl, Model& repModel, const CellParameters& parameters, bool verbose)
  {
    size_t numTrees = parameters.numTrees;
    double k        = parameters.k;
    double power    = parameters.power;
    double base     = parameters.base;
      
    if (verbose) ext_printf("    n.tree: %lu, k: %f, power: %f, base: %f\n", numTrees, k, power, base);
    
    repControl.numTrees = numTrees;
    
    double endNodeSd = (repControl.responseIsBinary ? 3.0 : 0.5) / (k * std::sqrt(static_cast<double>(repControl.numTrees)));
    static_cast<NormalPrior*>(repModel.muPrior)->precision = 1.0 / (endNodeSd * endNodeSd);
    
    static_cast<CGMPrior*>(repModel.treePrior)->power = power;
    static_cast<CGMPrior*>(repModel.treePrior)->base  = base;
    
    fit.setControl(repControl);
    fit.setModel(repModel);
  }
}

