#include "config.hpp"
#include "tree.hpp"

#include <cstring>
#include <cstdio>

#include <external/alloca.h>
#include <external/stats.h>
#include <external/linearAlgebra.h>

#include <dbarts/bartFit.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/scratch.hpp>
#include <dbarts/state.hpp>

namespace {
  using namespace dbarts;
  
  // multithread me!
  size_t* createObservationToNodeIndexMap(const BARTFit& fit, const Node& top,
                                          const double* xt, size_t numObservations)
  {
    if (numObservations == 0) return NULL;
    
    size_t* map = new size_t[numObservations];
        
    for (size_t i = 0; i < numObservations; ++i) {
      const Node* bottomNode = top.findBottomNode(fit, xt + i * fit.data.numPredictors);
      
      map[i] = bottomNode->enumerationIndex;
    }
    
    return map;
  }
}

// #include <external/io.h>
// #include <math.h>

namespace dbarts {
  void Tree::copyFrom(const BARTFit& fit, const Tree& other)
  {
    top.clear();
    
    std::memcpy(top.observationIndices, other.top.observationIndices, sizeof(size_t) * fit.data.numObservations);
    
    if (other.top.leftChild != NULL) {
      top.leftChild = new Node(top, fit.data.numPredictors, *other.top.leftChild);
      top.p.rightChild = new Node(top, fit.data.numPredictors, *other.top.p.rightChild);
      top.p.rule.copyFrom(other.top.p.rule);
    } else {
      top.m.average = other.top.m.average;
      top.m.numEffectiveObservations = other.top.m.numEffectiveObservations;
    }
    
    top.enumerationIndex = other.top.enumerationIndex;
    std::memcpy(top.variablesAvailableForSplit, other.top.variablesAvailableForSplit, sizeof(bool) * fit.data.numPredictors);
    
  }
  
  void Tree::setNodeAverages(const BARTFit& fit, size_t chainNum, const double* y) {
    NodeVector bottomNodes(getBottomNodes());
    
    size_t numBottomNodes = bottomNodes.size();
    
    for (size_t i = 0; i < numBottomNodes; ++i) {
      bottomNodes[i]->setAverage(fit, chainNum, y);
    }
  }
  
  void Tree::sampleAveragesAndSetFits(const BARTFit& fit, size_t chainNum, double sigma, double* trainingFits, double* testFits)
  {
    State& state(fit.state[chainNum]);
    
    NodeVector bottomNodes(top.getAndEnumerateBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    double* nodePosteriorPredictions = NULL;
    
    if (testFits != NULL) nodePosteriorPredictions = ext_stackAllocate(numBottomNodes, double);
    
    for (size_t i = 0; i < numBottomNodes; ++i) {
      const Node& bottomNode(*bottomNodes[i]);
      
      double posteriorPrediction = bottomNode.drawFromPosterior(state.rng, *fit.model.muPrior, sigma * sigma);
      bottomNode.setPredictions(trainingFits, posteriorPrediction);
      
      if (testFits != NULL) nodePosteriorPredictions[i] = posteriorPrediction;
    }
    
    if (testFits != NULL) {
      size_t* observationNodeMap = createObservationToNodeIndexMap(fit, top, fit.sharedScratch.xt_test, fit.data.numTestObservations);
      for (size_t i = 0; i < fit.data.numTestObservations; ++i) testFits[i] = nodePosteriorPredictions[observationNodeMap[i]];
      delete [] observationNodeMap;
      
      ext_stackFree(nodePosteriorPredictions);
    }
  }
  
  double* Tree::recoverAveragesFromFits(const BARTFit&, const double* treeFits)
  {
    NodeVector bottomNodes(top.getBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    double* result = new double[numBottomNodes];
    for (size_t i = 0; i < numBottomNodes; ++i) {
      if (bottomNodes[i]->isTop()) {
        result[i] = treeFits[0];
      } else if (bottomNodes[i]->getNumObservations() > 0) {
        result[i] = treeFits[bottomNodes[i]->observationIndices[0]];
      } else {
        result[i] = 0.0;
      }
    }

    return(result);
  }
  
  void Tree::setCurrentFitsFromAverages(const BARTFit& fit, const double* posteriorPredictions, double* trainingFits, double* testFits)
  {
    NodeVector bottomNodes(top.getAndEnumerateBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    if (trainingFits != NULL) {
      for (size_t i = 0; i < numBottomNodes; ++i) {
        const Node& bottomNode(*bottomNodes[i]);
        
        bottomNode.setPredictions(trainingFits, posteriorPredictions[i]);
      }
    }
    
    if (testFits != NULL) {
      size_t* observationNodeMap = createObservationToNodeIndexMap(fit, top, fit.sharedScratch.xt_test, fit.data.numTestObservations);
      for (size_t i = 0; i < fit.data.numTestObservations; ++i) testFits[i] = posteriorPredictions[observationNodeMap[i]];
      delete [] observationNodeMap;
    }
  }
  
  void Tree::setCurrentFitsFromAverages(const BARTFit& fit, const double* posteriorPredictions, const double* xt, size_t numObservations, double* fits)
  {
    top.enumerateBottomNodes();
    
    size_t* observationNodeMap = createObservationToNodeIndexMap(fit, top, xt, numObservations);
    for (size_t i = 0; i < numObservations; ++i) fits[i] = posteriorPredictions[observationNodeMap[i]];
    delete [] observationNodeMap;
  }
}

namespace {
  using namespace dbarts;
  void mapCutPoints(Node& n, const BARTFit& fit, const double* const* oldCutPoints, double* posteriorPredictions, int32_t* minIndices, int32_t* maxIndices, int32_t depth);
  void collapseEmptyNodes(Node& n, const BARTFit& fit, double* posteriorPredictions, int depth);
  void sampleFromPrior(const BARTFit& fit, ext_rng* rng, Node& n);
}

namespace dbarts {
  void Tree::mapOldCutPointsOntoNew(const BARTFit& fit, const double* const* oldCutPoints, double* posteriorPredictions)
  {
    // size_t origNumBottomNodes = top.getNumBottomNodes();
    
    int32_t* minIndices = new int32_t[fit.data.numPredictors];
    int32_t* maxIndices = new int32_t[fit.data.numPredictors];
    
    for (size_t i = 0; i < fit.data.numPredictors; ++i) {
      minIndices[i] = 0;
      maxIndices[i] = fit.sharedScratch.numCutsPerVariable[i];
    }
    
    mapCutPoints(top, fit, oldCutPoints, posteriorPredictions, minIndices, maxIndices, 2);
    
    delete [] maxIndices;
    delete [] minIndices;
   
    NodeVector bottomNodes(top.getBottomVector());
    size_t numBottomNodes = bottomNodes.size();
   
    for (size_t i = 0; i < numBottomNodes; ++i) {
      posteriorPredictions[i] = posteriorPredictions[bottomNodes[i]->enumerationIndex];
    }
  }
  
  void Tree::collapseEmptyNodes(const BARTFit& fit, double* posteriorPredictions)
  {
    // size_t origNumBottomNodes = top.getNumBottomNodes();
    
    top.enumerateBottomNodes();
    ::collapseEmptyNodes(top, fit, posteriorPredictions, 2);
    
    NodeVector bottomNodes(top.getBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    for (size_t i = 0; i < numBottomNodes; ++i) {
      posteriorPredictions[i] = posteriorPredictions[bottomNodes[i]->enumerationIndex];
    }
  }
  
  void Tree::countVariableUses(uint32_t* variableCounts) const {
    top.countVariableUses(variableCounts);
  }
  
  bool Tree::isValid() const {
    const NodeVector bottomNodes(top.getBottomVector());
    size_t numBottomNodes = bottomNodes.size();
    
    for (size_t j = 0; j < numBottomNodes; ++j) {
      if (bottomNodes[j]->getNumObservations() == 0) return false;
    }
    
    return true;
  }
  
  void Tree::sampleFromPrior(const BARTFit& fit, ext_rng* rng) {
    top.clear();
    ::sampleFromPrior(fit, rng, top);
  }
}

namespace {
  using namespace dbarts;
  
  // minIndex is inclusive, maxIndex is exclusive
  void mapCutPoints(Node& n, const BARTFit& fit, const double* const* oldCutPoints, double* posteriorPredictions, int32_t* minIndices, int32_t* maxIndices, int32_t depth)
  {
    if (n.isBottom() || n.p.rule.variableIndex == DBARTS_INVALID_RULE_VARIABLE) return;
    
    int32_t varIndex = n.p.rule.variableIndex;
    
    if (fit.data.variableTypes[varIndex] == ORDINAL) {
      int32_t minIndex = minIndices[varIndex];
      int32_t maxIndex = maxIndices[varIndex];
      
      double oldCut = oldCutPoints[varIndex][n.p.rule.splitIndex];
      const double* cutPoints_i = fit.sharedScratch.cutPoints[varIndex];
      
      
      if (minIndex > maxIndex - 1) {
        // no split can be made for this node, so we collapse it
        // since it is fundamentally invalid, we can't use a lot of information
        // from the nodes beneath it
        
        NodeVector bottomNodes(n.getBottomVector());
        size_t numBottomNodes = bottomNodes.size();
        double param = 0.0;
        for (size_t i = 0; i < numBottomNodes; ++i) param += posteriorPredictions[bottomNodes[i]->enumerationIndex];
        param /= static_cast<double>(numBottomNodes);
        
        size_t leftMostEnumerationIndex = bottomNodes[0]->enumerationIndex;
        delete n.getLeftChild();
        delete n.getRightChild();
        n.leftChild = NULL;
      
        posteriorPredictions[leftMostEnumerationIndex] = param;
        n.enumerationIndex = leftMostEnumerationIndex;
        return;
      } else {
        int32_t firstLessThan = n.p.rule.splitIndex < maxIndex ? n.p.rule.splitIndex : maxIndex - 1;
        // if it starts out below, move it above
        while (firstLessThan < maxIndex && cutPoints_i[firstLessThan] < oldCut) ++firstLessThan;
        // now nudge it back down
        if (firstLessThan < maxIndex) while (firstLessThan >= minIndex && cutPoints_i[firstLessThan] >= oldCut) --firstLessThan;
  
        int32_t newIndex;
        if (firstLessThan >= maxIndex - 1) newIndex = maxIndex - 1;
        else if (firstLessThan < minIndex) newIndex = minIndex;
        else if (cutPoints_i[firstLessThan + 1] == oldCut) newIndex = firstLessThan + 1;
        else if (oldCut - cutPoints_i[firstLessThan] < cutPoints_i[firstLessThan + 1] - oldCut) newIndex = firstLessThan;
        else newIndex = firstLessThan + 1;
        
        n.p.rule.splitIndex = newIndex;
      }
      
      maxIndices[varIndex] = n.p.rule.splitIndex;
      mapCutPoints(*n.leftChild, fit, oldCutPoints, posteriorPredictions, minIndices, maxIndices, depth + 1);
      maxIndices[varIndex] = maxIndex;
      
      minIndices[varIndex] = n.p.rule.splitIndex + 1;
      mapCutPoints(*n.p.rightChild, fit, oldCutPoints, posteriorPredictions, minIndices, maxIndices, depth + 1);
      minIndices[varIndex] = minIndex;
    }
  }
  
  void collapseEmptyNodes(Node& n, const BARTFit& fit, double* posteriorPredictions, int depth)
  {
    if (n.isBottom()) return; // only happens if is top and bottom
    
    if (n.getLeftChild()->getNumObservations() == 0 || n.getRightChild()->getNumObservations() == 0) {
      const NodeVector bottomNodes(n.getBottomVector());
      size_t numBottomNodes = bottomNodes.size();
      double* weights = ext_stackAllocate(numBottomNodes, double);
      double* params  = ext_stackAllocate(numBottomNodes, double);
      
      for (size_t i = 0; i < numBottomNodes; ++i) {
        Node& bottomNode(*bottomNodes[i]);
        weights[i] = fit.data.weights == NULL ? static_cast<double>(bottomNode.getNumObservations()) : ext_sumIndexedVectorElements(fit.data.weights, bottomNode.observationIndices, bottomNode.getNumObservations());
        params[i] = posteriorPredictions[bottomNodes[i]->enumerationIndex];
      }
      size_t leftMostEnumerationIndex = bottomNodes[0]->enumerationIndex;
      delete n.getLeftChild();
      delete n.getRightChild();
      n.leftChild = NULL;
      
      if (weights[0] == 0.0 && ext_vectorIsConstant(weights, numBottomNodes)) {
        posteriorPredictions[leftMostEnumerationIndex] = ext_computeMean(params, numBottomNodes);
      } else {
        posteriorPredictions[leftMostEnumerationIndex] = ext_computeWeightedMean(params, numBottomNodes, weights, NULL);
      }
      n.enumerationIndex = leftMostEnumerationIndex;
      
      ext_stackFree(params);
      ext_stackFree(weights);
    } else {
      if (!n.getLeftChild()->isBottom()) collapseEmptyNodes(*n.getLeftChild(), fit, posteriorPredictions, depth + 1);
      if (!n.getRightChild()->isBottom()) collapseEmptyNodes(*n.getRightChild(), fit, posteriorPredictions, depth + 1);
    }
  }
  
  void sampleFromPrior(const BARTFit& fit, ext_rng* rng, Node& n) {
    double parentPriorGrowthProbability = fit.model.treePrior->computeGrowthProbability(fit, n);
    if (parentPriorGrowthProbability <= 0.0 || ext_rng_simulateBernoulli(rng, parentPriorGrowthProbability) == 0) return;
    
    bool exhaustedLeftSplits, exhaustedRightSplits;
    Rule newRule = fit.model.treePrior->drawRuleAndVariable(fit, rng, n, &exhaustedLeftSplits, &exhaustedRightSplits);
    n.split(fit, newRule, exhaustedLeftSplits, exhaustedRightSplits);
    
    sampleFromPrior(fit, rng, *n.leftChild);
    sampleFromPrior(fit, rng, *n.p.rightChild);
  }
}
