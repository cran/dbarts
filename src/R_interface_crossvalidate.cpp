#include "config.hpp"
#include "R_interface_crossvalidate.hpp"

#include <cmath>   // sqrt, log
#include <cstddef> // size_t
#include <cstring> // strcmp

#include <misc/alloca.h>
#include <misc/linearAlgebra.h>
#include <misc/stats.h>
#include <misc/string.h>

#include <external/stats.h>

#include <R_ext/Random.h> // GetRNGstate, PutRNGState

#include <rc/bounds.h>
#include <rc/util.h>

#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>

#include "crossvalidate.hpp"

#include "R_interface_common.hpp"

#ifdef HAVE_LOG1P_IN_NAMESPACE_STD
using std::log1p;
#elif HAVE_LOG1P
#  include <math.h>
#else
double log1p(double x) {
  return std::log(1.0 + x);
}
#endif

using std::size_t;

namespace {
  using namespace dbarts;
  using namespace dbarts::xval;
  
  const char* const lossTypeNames[] = { "rmse", "log", "mcr", "custom" };
  typedef enum { RMSE, LOG, MCR, CUSTOM, INVALID } LossFunctorType;
  
  SEXP allocateResult(size_t numNTrees, size_t numKs, size_t numPowers, size_t numBases, size_t numReps, bool dropUnusedDims);
  LossFunctorDefinition* createLossFunctorDefinition(LossFunctorType lossType, SEXP lossTypeExpr, size_t numTestObservations,
                                                     size_t numSamples, SEXP scratch);
}

extern "C" {
  using namespace dbarts;
  using namespace dbarts::xval;


  SEXP xbart(SEXP controlExpr, SEXP modelExpr, SEXP dataExpr, SEXP methodExpr,
             SEXP testSampleSizeExpr, SEXP numRepsExpr, SEXP numBurnInExpr,
             SEXP lossTypeExpr, SEXP numThreadsExpr, 
             SEXP numTreesExpr, SEXP kExpr, SEXP powerExpr, SEXP baseExpr,
             SEXP dropExpr)
  {
    rc_assertIntConstraints(numTreesExpr, "num trees", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_VALUE | RC_GT, 0, RC_END);
    if (Rf_isReal(kExpr)) {
      rc_assertDoubleConstraints(kExpr, "k", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
    } else if (!Rf_isNull(kExpr)) {
      Rf_error("k must be numeric or NULL; hyperprior crossvalidation not supported at this time");
    }
       
    rc_assertDoubleConstraints(powerExpr, "power", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_END);
    rc_assertDoubleConstraints(baseExpr, "base", RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
    rc_assertIntConstraints(numBurnInExpr, "num burn",  RC_LENGTH | RC_GEQ, rc_asRLength(1), RC_LENGTH | RC_LEQ, rc_asRLength(3), RC_VALUE | RC_GEQ, 0, RC_END);
    
    Control control;
    Model model;
    Data data;
    
    SEXP classExpr = rc_getClass(controlExpr);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsControl") != 0) Rf_error("internal error: 'control' argument to dbarts_xbart not of class 'dbartsControl'");
    
    classExpr = rc_getClass(modelExpr);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsModel") != 0) Rf_error("internal error: 'model' argument to dbarts_xbart not of class 'dbartsModel'");
    
    classExpr = rc_getClass(dataExpr);
    if (std::strcmp(CHAR(STRING_ELT(classExpr, 0)), "dbartsData") != 0) Rf_error("internal error: 'data' argument to dbarts_xbart not of class 'dbartsData'");
    
    Method method;
    if (std::strcmp(CHAR(STRING_ELT(methodExpr, 0)), "k-fold") == 0) {
      method = K_FOLD;
    } else if (std::strcmp(CHAR(STRING_ELT(methodExpr, 0)), "random subsample") == 0) {
      method = RANDOM_SUBSAMPLE;
    } else {
      Rf_error("internal error: recognized method '%s'\n", CHAR(STRING_ELT(methodExpr, 0)));
    }
    
    // pull early so we don't allocate memory
    size_t numObservations = rc_getLength(Rf_getAttrib(dataExpr, Rf_install("y")));
    size_t numSamples      = static_cast<size_t>(INTEGER(Rf_getAttrib(controlExpr, Rf_install("n.samples")))[0]);
    
    sizetOrDouble testSampleSize;
    if (method == K_FOLD) 
      testSampleSize.n = rc_getInt(testSampleSizeExpr, "n.test", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 2, RC_VALUE | RC_LEQ, static_cast<int>(numObservations), RC_END);
    else
      testSampleSize.p = rc_getDouble(testSampleSizeExpr, "n.test", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0.0, RC_VALUE | RC_LT, 1.0, RC_END);
    
    size_t numReps = static_cast<size_t>(
      rc_getInt(numRepsExpr, "num reps", RC_LENGTH | RC_GEQ, rc_asRLength(1),
                                         RC_VALUE | RC_GT, 0, RC_END));
    
    int numThreadsInt = rc_getInt(numThreadsExpr, "num threads", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_VALUE | RC_GT, 0, RC_NA | RC_YES, RC_END);
    size_t numThreads = numThreadsInt != NA_INTEGER ? static_cast<size_t>(numThreadsInt) : 1;

    size_t numInitialBurnIn      = static_cast<size_t>(INTEGER(numBurnInExpr)[0]);
    size_t numContextShiftBurnIn = rc_getLength(numBurnInExpr) >= 2 ? static_cast<size_t>(INTEGER(numBurnInExpr)[1]) : ((3 * numInitialBurnIn) / 4);
    size_t numRepBurnIn          = rc_getLength(numBurnInExpr) == 3 ? static_cast<size_t>(INTEGER(numBurnInExpr)[2]) : numInitialBurnIn / 4;
    
    bool dropUnusedDims = rc_getBool(dropExpr, "drop", RC_LENGTH | RC_EQ, rc_asRLength(1), RC_END);
    
    size_t maxNumTestObservations;
    if (method == K_FOLD) {
      maxNumTestObservations = numObservations / testSampleSize.n;
      if (numObservations % maxNumTestObservations != 0) ++maxNumTestObservations;
    } else {
      maxNumTestObservations = numObservations -
        static_cast<size_t>(std::floor(static_cast<double>(numObservations) * testSampleSize.p + 0.5));
    }
    
    LossFunctorType lossType = INVALID;
    if (Rf_isString(lossTypeExpr)) {
      if (rc_getLength(lossTypeExpr) != 1) Rf_error("length of lossType for strings must be 1");
      const char* lossTypeName = CHAR(STRING_ELT(lossTypeExpr, 0));
      
      size_t lossTypeNumber;
      int errorCode = misc_str_matchInArray(lossTypeName, lossTypeNames, static_cast<size_t>(CUSTOM - RMSE), &lossTypeNumber);
      if (errorCode != 0) Rf_error("error matching string: %s", std::strerror(errorCode));
      if (lossTypeNumber == MISC_STR_NO_MATCH) Rf_error("unsupported result type: '%s'", lossTypeName);
      
      lossType = static_cast<LossFunctorType>(lossTypeNumber);
    } else if (Rf_isVectorList(lossTypeExpr)) {
      if (rc_getLength(lossTypeExpr) != 2) Rf_error("length of lossType for functions must be 2");
      
      if (!Rf_isFunction(VECTOR_ELT(lossTypeExpr, 0))) Rf_error("first element of list for function lossType must be a closure");
      if (!Rf_isEnvironment(VECTOR_ELT(lossTypeExpr, 1))) Rf_error("second element of list for function lossType must be an environment");
      
      lossType = CUSTOM;
    } else {
      Rf_error("lossType must be a character string or list(closure, env)");
    }
    
    int protectCount = 0;
    
    SEXP scratch = R_NilValue;
    if (lossType == CUSTOM) {
      R_xlen_t scratchLength = rc_asRLength(numThreads * 3 * (method == K_FOLD ? 2 : 1));
      scratch = PROTECT(rc_newList(scratchLength));
      ++protectCount;
      for (R_xlen_t i = 0; i < scratchLength; ++i)
        SET_VECTOR_ELT(scratch, i, R_NilValue);
    }
    
    initializeControlFromExpression(control, controlExpr);
    if (control.defaultNumSamples == 0) {
      if (protectCount > 0) UNPROTECT(protectCount);
      Rf_error("xbart called with 0 posterior samples");
    }
    
    if (lossType == LOG && control.responseIsBinary == false) lossType = RMSE;
    
    LossFunctorDefinition* lossFunctionDef =
      createLossFunctorDefinition(lossType, lossTypeExpr, maxNumTestObservations, numSamples, scratch);
    
    initializeModelFromExpression(model, modelExpr, control);
    initializeDataFromExpression(data, dataExpr);
    
    if (data.numObservations == 0) {
      invalidateData(data);
      invalidateModel(model);
      delete lossFunctionDef;
      
      if (protectCount > 0) UNPROTECT(protectCount);
      Rf_error("xbart called on empty data set");
    }
    
    size_t numNTrees = rc_getLength(numTreesExpr);
    size_t numKs     = Rf_isReal(kExpr) ? rc_getLength(kExpr) : 1;
    size_t numPowers = rc_getLength(powerExpr);
    size_t numBases  = rc_getLength(baseExpr);
    
    int* nTreesInt  = INTEGER(numTreesExpr);
    size_t* nTrees = misc_stackAllocate(numNTrees, size_t);
    for (size_t i = 0; i < numNTrees; ++i) nTrees[i] = static_cast<size_t>(nTreesInt[i]);
    
    double* k     = Rf_isReal(kExpr) ? REAL(kExpr) : NULL;
    double* power = REAL(powerExpr);
    double* base  = REAL(baseExpr);
    
    SEXP result = PROTECT(allocateResult(numNTrees, numKs, numPowers, numBases, numReps, dropUnusedDims));
    ++protectCount;
    
    GetRNGstate();
    
    crossvalidate(control, model, data, method,
                  testSampleSize, numReps,
                  numInitialBurnIn, numContextShiftBurnIn, numRepBurnIn,
                  *lossFunctionDef, numThreads,
                  nTrees, numNTrees, k, numKs, power, numPowers, base, numBases,
                  REAL(result));
    
    PutRNGstate();
    
    delete lossFunctionDef;
    
    misc_stackFree(nTrees);
    
    invalidateData(data);
    invalidateModel(model);

    UNPROTECT(protectCount);
    
    return result;
  }
}

namespace {
  using namespace dbarts;
  using namespace dbarts::xval;
  
  SEXP allocateResult(size_t numNTrees, size_t numKs, size_t numPowers, size_t numBases, size_t numReps, bool dropUnusedDims)
  {
    SEXP result = PROTECT(Rf_allocVector(REALSXP, numNTrees * numKs * numPowers * numBases * numReps));
    if (dropUnusedDims) {
      size_t numDims = 1 + (numNTrees > 1 ? 1 : 0) + (numKs > 1 ? 1 : 0) + (numPowers > 1 ? 1 : 0) + (numBases > 1 ? 1 : 0);
      if (numDims > 1) {
        SEXP dimsExpr = PROTECT(Rf_allocVector(INTSXP, numDims));
        int* dims = INTEGER(dimsExpr);
        dims[0] = static_cast<int>(numReps);
        numDims = 1;
        if (numNTrees > 1) dims[numDims++] = static_cast<int>(numNTrees);
        if (numKs > 1)     dims[numDims++] = static_cast<int>(numKs);
        if (numPowers > 1) dims[numDims++] = static_cast<int>(numPowers);
        if (numBases > 1)  dims[numDims]   = static_cast<int>(numBases);
        
        R_do_slot_assign(result, R_DimSymbol, dimsExpr);
        UNPROTECT(1);
      }
    } else {
      rc_setDims(result, static_cast<int>(numReps), static_cast<int>(numNTrees), static_cast<int>(numKs), static_cast<int>(numPowers),
                 static_cast<int>(numBases), -1);
    }
    UNPROTECT(1);
    return result;
  }
  
  struct LogLossFunctor : LossFunctor {
    double* scratch;
  };
  
  LossFunctor* createLogLoss(const LossFunctorDefinition& def, Method, std::size_t numTestObservations, std::size_t numSamples)
  {
    (void) def; (void) numSamples; (void) numTestObservations;
    
    LogLossFunctor* result = new LogLossFunctor;
    result->scratch = new double[numSamples];
    return result;
  }
  
  void deleteLogLoss(LossFunctor* instance)
  {
    delete [] static_cast<LogLossFunctor*>(instance)->scratch;
    delete static_cast<LogLossFunctor*>(instance);
  }
  
  void calculateLogLoss(LossFunctor& restrict v_instance,
                        const double* restrict y_test, size_t numTestObservations, const double* restrict testSamples, size_t numSamples,
                        double* restrict results)
  {
    LogLossFunctor& restrict instance(*static_cast<LogLossFunctor* restrict>(&v_instance));
    
    double* restrict probabilities = instance.scratch;
    results[0] = 0.0;
    for (size_t i = 0; i < numTestObservations; ++i) {
      for (size_t j = 0; j < numSamples; ++j) {
        probabilities[j] = ext_cumulativeProbabilityOfNormal(testSamples[i + j * numTestObservations], 0.0, 1.0);
      }
      double y_test_hat = misc_computeMean(probabilities, numSamples);
      
      results[0] +=  -y_test[i] * std::log(y_test_hat) - (1.0 - y_test[i]) * log1p(-y_test_hat);
    }
    results[0] /= static_cast<double>(numTestObservations);
  }
  
  struct RMSELossFunctor : LossFunctor {
    double* scratch;
  };
  
  LossFunctor* createRMSELoss(const LossFunctorDefinition& def, Method, std::size_t numTestObservations, std::size_t numSamples)
  {
    (void) def; (void) numSamples;
    
    RMSELossFunctor* result = new RMSELossFunctor;
    result->scratch = new double[numTestObservations];
    return result;
  }
  
  void deleteRMSELoss(LossFunctor* instance)
  {
    delete [] static_cast<RMSELossFunctor*>(instance)->scratch;
    delete static_cast<RMSELossFunctor*>(instance);
  }
  
  void calculateRMSELoss(LossFunctor& restrict v_instance,
                        const double* restrict y_test, size_t numTestObservations, const double* restrict testSamples, size_t numSamples,
                        double* restrict results)
  {
    RMSELossFunctor& restrict instance(*static_cast<RMSELossFunctor* restrict>(&v_instance));
    
    double* restrict y_test_hat = instance.scratch;
    for (size_t i = 0; i < numTestObservations; ++i) {
      y_test_hat[i] = 0.0;
      for (size_t j = 0; j < numSamples; ++j) y_test_hat[i] += testSamples[i + j * numTestObservations];
      y_test_hat[i] /= static_cast<double>(numSamples);
    }
    
    results[0] = std::sqrt(misc_computeSumOfSquaredResiduals(y_test, numTestObservations, y_test_hat) / static_cast<double>(numTestObservations));
  }

  
  struct MCRLossFunctor : LossFunctor {
    double* scratch;
  };
  
  LossFunctor* createMCRLoss(const LossFunctorDefinition& def, Method, std::size_t numTestObservations, std::size_t numSamples)
  {
    (void) def; (void) numTestObservations;
    
    MCRLossFunctor* result = new MCRLossFunctor;
    result->scratch = new double[numSamples];
    return result;
  }
  
  void deleteMCRLoss(LossFunctor* instance)
  {
    delete [] static_cast<MCRLossFunctor*>(instance)->scratch;
    delete static_cast<MCRLossFunctor*>(instance);
  }

  void calculateMCRLoss(LossFunctor& restrict v_instance,
                        const double* restrict y_test, size_t numTestObservations, const double* restrict testSamples, size_t numSamples,
                        double* restrict results)
  {
    MCRLossFunctor& restrict instance(*static_cast<MCRLossFunctor* restrict>(&v_instance));
    
    double* restrict probabilities = instance.scratch;
    size_t fp = 0, fn = 0;
    for (size_t i = 0; i < numTestObservations; ++i) {
      for (size_t j = 0; j < numSamples; ++j) {
        probabilities[j] = ext_cumulativeProbabilityOfNormal(testSamples[i + j * numTestObservations], 0.0, 1.0);
      }
      double y_test_hat = misc_computeMean(probabilities, numSamples) > 0.5 ? 1.0 : 0.0;
      
      if (y_test[i] != y_test_hat) {
        if (y_test[i] == 1.0) ++fn; else ++fp;
      }
    }
    results[0] = static_cast<double>(fp + fn) / static_cast<double>(numTestObservations);
  }
  
  /*
   * Custom loss is created from a function and the environment to evaluate that function in.
   * In order to call the function, we allocate and store R vectors for y_test, samples 
   * of y_test, and the closure of the function applied to those values. They are released
   * to the garbage collector when the loss function is deleted.
   *
   * Since K-Fold crossvalidation can have blocks with unequal test sample sizes, we
   * create extra storage and a closure to handle.
   */
  struct CustomLossFunctorDefinition : LossFunctorDefinition {
    SEXP function;
    SEXP environment;
    SEXP scratch;
        
    ~CustomLossFunctorDefinition() { }
  };
  
  struct CustomLossFunctor : LossFunctor {
    double* y_test;
    double* testSamples;
    size_t maxNumTestObservations;
    
    double* y_testNM1;
    double* testSamplesNM1;
    
    SEXP closure;
    SEXP closureNM1;
    SEXP environment;
  };
  
  LossFunctor* createCustomLoss(const LossFunctorDefinition& v_def, Method method, std::size_t numTestObservations, std::size_t numSamples)
  {
    CustomLossFunctorDefinition& def(*const_cast<CustomLossFunctorDefinition*>(static_cast<const CustomLossFunctorDefinition*>(&v_def)));
    CustomLossFunctor* result = new CustomLossFunctor;
    
    size_t assignmentIndex;
    size_t scratchLength = rc_getLength(def.scratch);
    for (assignmentIndex = 0; assignmentIndex < scratchLength && VECTOR_ELT(def.scratch, assignmentIndex) != R_NilValue; ++assignmentIndex) { /* */ }
    
    SEXP y_testExpr      = PROTECT(Rf_allocVector(REALSXP, numTestObservations));
    SEXP testSamplesExpr = PROTECT(Rf_allocVector(REALSXP, numTestObservations * numSamples));
    rc_setDims(testSamplesExpr, static_cast<int>(numTestObservations), static_cast<int>(numSamples), -1);
            
    result->y_test      = REAL(y_testExpr);
    result->testSamples = REAL(testSamplesExpr);
    result->maxNumTestObservations = numTestObservations;
    
    result->closure     = PROTECT(Rf_lang3(def.function, y_testExpr, testSamplesExpr));
    result->environment = def.environment;
    
    SET_VECTOR_ELT(def.scratch, assignmentIndex++, y_testExpr);
    SET_VECTOR_ELT(def.scratch, assignmentIndex++, testSamplesExpr);
    SET_VECTOR_ELT(def.scratch, assignmentIndex++, result->closure);
    
    UNPROTECT(3);
    
    if (method == K_FOLD) {
      y_testExpr      = PROTECT(Rf_allocVector(REALSXP, numTestObservations - 1));
      testSamplesExpr = PROTECT(Rf_allocVector(REALSXP, (numTestObservations - 1) * numSamples));
      rc_setDims(testSamplesExpr, static_cast<int>(numTestObservations - 1), static_cast<int>(numSamples), -1);
      
      result->y_testNM1      = REAL(y_testExpr);
      result->testSamplesNM1 = REAL(testSamplesExpr);
            
      result->closureNM1 = PROTECT(Rf_lang3(def.function, y_testExpr, testSamplesExpr));
      
      SET_VECTOR_ELT(def.scratch, assignmentIndex++, y_testExpr);
      SET_VECTOR_ELT(def.scratch, assignmentIndex++, testSamplesExpr);
      SET_VECTOR_ELT(def.scratch, assignmentIndex, result->closureNM1);
      
      UNPROTECT(3);
    }
    
    return result;
  }
  
  void deleteCustomLoss(LossFunctor* v_instance)
  {
    delete static_cast<CustomLossFunctor*>(v_instance);
  }
  
  void calculateCustomLoss(LossFunctor& restrict v_instance,
                           const double* restrict, size_t numTestObservations, const double* restrict, size_t numSamples,
                           double* restrict results)
  {
    CustomLossFunctor& restrict instance(*static_cast<CustomLossFunctor* restrict>(&v_instance));
    
    SEXP customResult;
    if (numTestObservations == instance.maxNumTestObservations) {
      customResult = Rf_eval(instance.closure, instance.environment);
    } else {
      std::memcpy(instance.y_testNM1, const_cast<const double*>(instance.y_test), numTestObservations * sizeof(double));
      std::memcpy(instance.testSamplesNM1, const_cast<const double*>(instance.testSamples), numTestObservations * numSamples * sizeof(double));
      customResult = Rf_eval(instance.closureNM1, instance.environment);
    }
    
    std::memcpy(results, const_cast<const double*>(REAL(customResult)), rc_getLength(customResult) * sizeof(double));
  }
  
  LossFunctorDefinition* createLossFunctorDefinition(LossFunctorType type, SEXP lossTypeExpr, size_t numTestObservations, size_t numSamples,
                                                     SEXP scratch)
  {
    LossFunctorDefinition* result = NULL;
    
    switch(type) {
      case RMSE:
      result = new LossFunctorDefinition;
      result->y_testOffset = -1;
      result->testSamplesOffset = -1;
      result->numResults = 1;
      result->displayString = lossTypeNames[RMSE];
      result->requiresMutex = false;
      result->calculateLoss = &calculateRMSELoss;
      result->createFunctor = &createRMSELoss;
      result->deleteFunctor = &deleteRMSELoss;
      break;
      case LOG:
      result = new LossFunctorDefinition;
      result->y_testOffset = -1;
      result->testSamplesOffset = -1;
      result->numResults = 1;
      result->displayString = lossTypeNames[LOG];
      result->requiresMutex = false;
      result->calculateLoss = &calculateLogLoss;
      result->createFunctor = &createLogLoss;
      result->deleteFunctor = &deleteLogLoss;
      break;
      case MCR:
      result = new LossFunctorDefinition;
      result->y_testOffset = -1;
      result->testSamplesOffset = -1;
      result->numResults = 1;
      result->displayString = lossTypeNames[MCR];
      result->requiresMutex = false;
      result->calculateLoss = &calculateMCRLoss;
      result->createFunctor = &createMCRLoss;
      result->deleteFunctor = &deleteMCRLoss;
      break;
      case CUSTOM:
      {
        SEXP function    = VECTOR_ELT(lossTypeExpr, 0);
        SEXP environment = VECTOR_ELT(lossTypeExpr, 1);
        
        CustomLossFunctorDefinition* c_result = new CustomLossFunctorDefinition;
        result = c_result;
        result->y_testOffset = 0;
        result->testSamplesOffset = sizeof(double*);
        
        SEXP tempY_test      = PROTECT(Rf_allocVector(REALSXP, numTestObservations));
        SEXP tempTestSamples = PROTECT(Rf_allocVector(REALSXP, numTestObservations * numSamples));
        rc_setDims(tempTestSamples, static_cast<int>(numTestObservations), static_cast<int>(numSamples), -1);
        // set some values to not all be the same so that they have an SD != 0
        REAL(tempY_test)[0] = 1.0;
        misc_setVectorToConstant(REAL(tempY_test) + 1, numTestObservations - 1, 0.0);
        REAL(tempTestSamples)[0] = 0.0;
        misc_setVectorToConstant(REAL(tempTestSamples) + 1, numTestObservations - 1, 1.0);
        misc_setVectorToConstant(REAL(tempTestSamples) + numTestObservations, numTestObservations * (numSamples - 1), 0.0);
        
        SEXP tempClosure = PROTECT(Rf_lang3(function, tempY_test, tempTestSamples));
        
        result->numResults = rc_getLength(Rf_eval(tempClosure, environment));
        UNPROTECT(3);
        
        result->displayString = lossTypeNames[CUSTOM];
        result->requiresMutex = true;
        result->calculateLoss = &calculateCustomLoss;
        result->createFunctor = &createCustomLoss;
        result->deleteFunctor = &deleteCustomLoss;
        
        c_result->function = function;
        c_result->environment = environment;
        c_result->scratch = scratch;
      }
      break;
      case INVALID:
      Rf_error("internal error: invalid type enumeration");
      break;
    }
    
    return result;
  }
}

