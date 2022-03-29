## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----simulateData-------------------------------------------------------------
f <- function(x)
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]

set.seed(99)
sigma <- 1.0
n     <- 100

x <- matrix(runif(n * 10), n, 10)
y <- rnorm(n, f(x), sigma)

data <- data.frame(x, y)

## ----fitModel-----------------------------------------------------------------
library(dbarts, quietly = TRUE)

bartFit <- bart(
    y ~ ., data,
    ndpost = 4,   # number of posterior samples
    nskip = 1000, # number of "warmup" samples to discard
    nchain = 2,   # number of independent, parallel chains
    nthread = 1,  # units of parallel execution
    ntree = 3,    # number of trees per chain
    seed = 2,     # chosen to generate a deep tree
    keeptrees = TRUE,
    verbose = FALSE)

## ----extractTrees-------------------------------------------------------------
trees <- extract(bartFit, "trees")

## ----printFlattenedTrees------------------------------------------------------
print(head(trees, n = 10))

## ----rebuildTrees-------------------------------------------------------------
# Turns a flatted tree data frame into a list of lists, or a "natural" tree
# structure.
rebuildTree <- function(tree, object) {
    # Define a worker function that will be recursively called on every node.
    rebuildTreeRecurse <- function(tree) {
        node <- list(
            value = tree$value[1],
            n     = tree$n[1]
        )
        # Check node if is a leaf, and if so return early.
        if (tree$var[1] == -1) {
            node$n_nodes <- 1
            return(node)
        }
        
        node$var <- variableNames[tree$var[1]]
        
        # By removing the current row, we can recurse down the left branch.
        headOfLeftBranch <- tree[-1,]
        left <- rebuildTreeRecurse(headOfLeftBranch)
        n_nodes.left <- left$n_nodes
        left$n_nodes <- NULL
        node$left <- left
        
        # The right branch is obtained by advancing past the left nodes.
        headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
        right <- rebuildTreeRecurse(headOfRightBranch)
        n_nodes.right <- right$n_nodes
        right$n_nodes <- NULL
        node$right <- right
        
        node$n_nodes <- 1L + n_nodes.left + n_nodes.right
        
        return(node)
    }
    variableNames <- colnames(object$fit$data@x)
    
    result <- rebuildTreeRecurse(tree)
    result$n_nodes <- NULL
    return(result)
}

treeOfInterest <- subset(trees, chain == 1 & sample == 3 & tree == 1)
print(rebuildTree(treeOfInterest, bartFit))

## ----rebuildAllTrees----------------------------------------------------------
allTrees <- by(
    data    = trees,
    INDICES = trees[,c("chain", "sample", "tree")],
    FUN     = rebuildTree, 
    object  = bartFit)

# One way to index the result of this:
#    allTrees[chain = "1", sample = "2", tree = "3"]

## ----plotTree-----------------------------------------------------------------
bartFit$fit$plotTree(chainNum = 1, sampleNum = 3, treeNum = 1)

## ----getPredictionsForTree----------------------------------------------------
getPredictionsForTree <- function(tree, x) {
    predictions <- rep(NA_real_, nrow(x))
    
    getPredictionsForTreeRecursive <- function(tree, indices) {
        if (tree$var[1] == -1) {
            # Assigns in the calling environment by using <<-
            predictions[indices] <<- tree$value[1]
            return(1)
        }

        goesLeft <- x[indices, tree$var[1]] <= tree$value[1]
        headOfLeftBranch <- tree[-1,]
        n_nodes.left <- getPredictionsForTreeRecursive(
            headOfLeftBranch, indices[goesLeft])
        
        headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
        n_nodes.right <- getPredictionsForTreeRecursive(
            headOfRightBranch, indices[!goesLeft])
        
        return(1 + n_nodes.left + n_nodes.right)
    }

    getPredictionsForTreeRecursive(tree, seq_len(nrow(x)))

    return(predictions)
}

getPredictionsForTree(treeOfInterest, bartFit$fit$data@x[1:5,])

## ----mapOverNodes-------------------------------------------------------------
mapOverNodes <- function(tree, f, ...) {
    mapOverNodesRecurse <- function(tree, depth, f, ...) {
        node <- list(
            value = tree$value[1],
            n = tree$n[1],
            depth = depth
        )
        if (tree$var[1] == -1) {
            node$n_nodes <- 1
            node$f.x <- f(node, ...)
            return(node)
        }
        node$var <- tree$var[1]
        node$f.x <- f(node, ...)
        
        headOfLeftBranch <- tree[-1,]
        left <- mapOverNodesRecurse(headOfLeftBranch, depth + 1, f, ...)
        n_nodes.left <- left$n_nodes
        left$n_nodes <- NULL
        node$left <- left

        
        headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
        right <- mapOverNodesRecurse(headOfRightBranch, depth + 1, f, ...)
        n_nodes.right <- right$n_nodes
        right$n_nodes <- NULL
        node$right <- right
        
        node$n_nodes <- 1 + n_nodes.left + n_nodes.right
        return(node)
    }
    result <- mapOverNodesRecurse(tree, 1, f, ...)
    result$n_nodes <- NULL
    return(result)
}

## ----observeInteractions------------------------------------------------------
observeInteractions <- function(node, ...) {
    if (is.null(node$var)) return(NULL)

    interactionData <- list(...)$interactionData
    # Make the current node visibile inside the environment.
    interactionData$node <- node
    with(interactionData, {
        if (node$depth <= currentDepth) {
            # If true, we have backtracked to go down the right branch, so we
            # remove the variables from the left branch.
            currentVariables <- currentVariables[seq_len(node$depth - 1)]
        }
        if (length(interactionData$currentVariables) > 0) {
            # This is a brute-force way of updating the following indices,
            # relying on the column-major storage order that R uses:
            #     hasInteraction[currentVariables,,drop = FALSE][,node$var]
            updateIndices <- currentVariables +
                (node$var - 1) * nrow(hasInteraction)
            hasInteraction[updateIndices] <- TRUE
        }
        currentVariables <- c(currentVariables, node$var)
        currentDepth <- node$depth
    })
    rm("node", envir = interactionData)
    
    # Since the function is used for its side effects, there isn't a return
    # value.
    return(NULL)
}

numVariables  <- ncol(bartFit$fit$data@x)
variableNames <- colnames(bartFit$fit$data@x)

# Define this as an environment as they are mutable
interactionData <- list2env(list(
    currentDepth = 0,
    currentVariables = integer(),
    hasInteraction = matrix(
        data = FALSE,
        ncol = numVariables, nrow = numVariables,
        dimnames = list(ancestor = variableNames, descendant = variableNames)
    )
))

invisible(mapOverNodes(
    treeOfInterest,
    observeInteractions,
    interactionData = interactionData
))

## ----printObserveInteractionResults-------------------------------------------
print(interactionData$hasInteraction)

