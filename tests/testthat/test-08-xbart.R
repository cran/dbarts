context("xbart")

source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

test_that("random subsample runs correctly with valid inputs", {
  x <- testData$x
  y <- testData$y
  
  n.reps  <- 3L
  n.trees <- c(5L, 7L)
  k       <- c(1, 2, 4)
  power   <- c(1.5, 2)
  base    <- c(0.75, 0.8, 0.95)

  xval <- xbart(x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "random subsample",
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base,
                n.threads = 2L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
})

test_that("k-fold runs correctly with valid inputs", {
  x <- testData$x
  y <- testData$y
  
  n.reps  <- 3L
  n.trees <- c(5L, 7L)
  k       <- c(1, 2, 4)
  power   <- c(1.5, 2)
  base    <- c(0.75, 0.8, 0.95)

  xval <- xbart(x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base,
                n.threads = 2L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
})

test_that("k-fold runs correctly with one input", {
  x <- testData$x
  y <- testData$y
  
  xval <- xbart(x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L),
                n.reps = 3, n.test = 5,
                k = 2,
                n.threads = 2L)
  
  expect_equal(length(xval), 3)
})

test_that("runs correctly with weighted input", {
  x <- testData$x
  y <- testData$y
  weights <- rep(1, length(y))
  
  xval <- xbart(x, y, weights = weights,
                n.samples = 6L, n.burn = c(5L, 3L, 1L),
                n.reps = 3, n.test = 5,
                k = 2,
                n.threads = 2L)
  
  expect_equal(length(xval), 3)
})


test_that("k-fold and random subsample are reproducible, roughly similar", {
  x <- testData$x
  y <- testData$y
  
  k <- c(4, 8)
  
  set.seed(0)
  xval.kf <- xbart(x, y, method = "k-fold",
                   n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
                   k = k,
                   n.threads = 1L)

  xval.rs <- xbart(x, y, method = "random subsample", 
                   n.reps = 20L, n.samples = 20L, n.burn = c(10L, 5L, 1L),
                   k = k, 
                   n.threads = 1L)
  
  res.rs <- apply(xval.rs, 2, mean)
  res.kf <- apply(xval.kf, 2, mean)
  expect_equal(unname(res.rs), c(2.35131034003841, 4.57299444101639))
  expect_equal(unname(res.kf), c(2.30094811299725, 4.54475202324197))
  
  expect_true(all(abs(res.rs - res.kf) < .1))
})

test_that("works with fixed seed", {
  x <- testData$x
  y <- testData$y
  
  k <- c(4, 8)
  
  xval.1 <- xbart(x, y, method = "k-fold",
                  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
                  k = k,
                  n.threads = 1L, seed = 0)

  xval.2 <- xbart(x, y, method = "k-fold",
                  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
                  k = k,
                  n.threads = 1L, seed = 0)
  
  expect_true(all(!is.na(xval.1)))
  expect_equal(dim(xval.1), c(4L, length(k)))
  expect_equal(xval.1, xval.2)

  xval.3 <- xbart(x, y, method = "k-fold",
                  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
                  k = k,
                  n.threads = 2L, seed = 0)

  xval.4 <- xbart(x, y, method = "k-fold",
                  n.reps = 4L, n.samples = 20L, n.burn = c(10L, 5L, 1L), n.test = 5,
                  k = k,
                  n.threads = 2L, seed = 0)
  
  expect_true(all(!is.na(xval.3)))
  expect_equal(dim(xval.3), c(4L, length(k)))
  expect_equal(xval.3, xval.4)

  expect_true(any(xval.1 != xval.3))
})

test_that("works with non-standard models", {
  x <- testData$x
  y <- testData$y
  
  k <- c(4, 8)
  expect_silent(xbart(x, y, method = "k-fold",
                      n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
                      k = k, n.threads = 1L, resid.prior = chisq(2.5, 0.9)))
  
  expect_silent(xbart(x, y, method = "k-fold",
                      n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
                      k = k, n.threads = 1L, resid.prior = fixed(2)))
  
  n.trees <- c(5, 10)
  expect_silent(xbart(x, y, method = "k-fold",
                      n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
                      n.trees = n.trees, n.threads = 1L))
})

test_that("works with custom loss", {
  x <- testData$x
  y <- testData$y
  
  n.reps  <- 3L
  n.trees <- c(5L, 7L)
  k       <- c(1, 2, 4)
  power   <- c(1.5, 2)
  base    <- c(0.75, 0.8, 0.95)
  
  mad <- function(y.train, y.train.hat, weights) 
    mean(abs(y.train - apply(y.train.hat, 1L, mean)))

  xval <- xbart(x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base, loss = mad,
                n.threads = 1L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
  
  
  xval <- xbart(x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base, loss = mad,
                n.threads = 2L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
})

test_that("fails with invalid inputs", {
  x <- testData$x
  y <- testData$y
  
  expect_error(xbart(y ~ x, method = "not-a-method"))
  expect_error(xbart(y ~ x, method = NULL))
  expect_error(xbart(y ~ x, method = NA_character_))
  
  expect_error(xbart(y ~ x, n.samples = 0L))
  expect_error(xbart(y ~ x, n.samples = "not-a-integer"))
  expect_error(xbart(y ~ x, n.samples = NULL))
  expect_error(xbart(y ~ x, n.samples = NA_integer_))
  
  expect_error(xbart(y ~ x, method = "k-fold", n.test = 1))
  expect_error(xbart(y ~ x, method = "k-fold", n.test = length(testData$y) + 1))
  expect_error(xbart(y ~ x, method = "random subsample", n.test = 0))
  expect_error(xbart(y ~ x, method = "random subsample", n.test = length(testData$y) + 1))
  expect_error(xbart(y ~ x, n.test = "not-a-numeric"))
  expect_error(xbart(y ~ x, n.test = NULL))
  expect_error(xbart(y ~ x, n.test = NA_real_))
  
  expect_error(xbart(y ~ x, n.reps = 0L))
  expect_error(xbart(y ~ x, n.reps = "not-a-integer"))
  expect_error(xbart(y ~ x, n.reps = NULL))
  expect_error(xbart(y ~ x, n.reps = NA_integer_))
  
  expect_error(xbart(y ~ x, n.burn = c(200L, -1L, 50L)))
  expect_error(xbart(y ~ x, n.burn = "not-a-integer"))
  expect_error(xbart(y ~ x, n.burn = NULL))
  expect_error(xbart(y ~ x, n.burn = NA_integer_))
  
  expect_error(xbart(y ~ x, loss = "unknown-loss"))
  expect_error(xbart(y ~ x, loss = 2))
  expect_error(xbart(y ~ x, loss = function(x) x))
  expect_error(xbart(y ~ x, loss = list(function(x, y) NULL, "not-an-environment")))
  
  expect_error(xbart(y ~ x, n.threads = 0L))
  expect_error(xbart(y ~ x, n.threads = "not-a-integer"))
  expect_error(xbart(y ~ x, n.threads = NULL))
  
  expect_error(xbart(y ~ x, n.trees = 0L))
  expect_error(xbart(y ~ x, n.trees = "not-a-integer"))
  expect_error(xbart(y ~ x, n.trees = NULL))
  expect_error(xbart(y ~ x, n.trees = NA_integer_))
  
  expect_error(xbart(y ~ x, k = c(-0.5, 1)))
  expect_error(xbart(y ~ x, k = "not-a-numeric"))
  expect_error(xbart(y ~ x, k = NA_real_))
  
  expect_error(xbart(y ~ x, power = c(0, 0.5)))
  expect_error(xbart(y ~ x, power = "not-a-numeric"))
  expect_error(xbart(y ~ x, power = NULL))
  expect_error(xbart(y ~ x, power = NA_real_))
  
  expect_error(xbart(y ~ x, base = c(0.5, 1)))
  expect_error(xbart(y ~ x, base = "not-a-numeric"))
  expect_error(xbart(y ~ x, base = NULL))
  expect_error(xbart(y ~ x, base = NA_real_))
  
  expect_error(xbart(y ~ x, verbose = "not-a-logical"))
  expect_error(xbart(y ~ x, verbose = NULL))
  expect_error(xbart(y ~ x, verbose = NA))
  
  expect_error(xbart(y ~ x, resid.prior = NULL))
  
  expect_error(xbart(y ~ x, sigma = -1))
  expect_error(xbart(y ~ x, sigma = "not-a-numeric"))
})

test_that("k-fold subdivides data correctly when data do not divide evenly by k", {
  x <- testData$x[1:24,]
  y <- testData$y[1:24]
  
  k <- c(2, 4)

  xval <- xbart(x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
                n.reps = 3L,
                k = k,
                n.threads = 1L)
  
  expect_is(xval, "array")
})

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

test_that("runs with binary data and k hyperprior", {
  x <- testData$X
  z <- testData$Z
  
  n.reps  <- 3L
  power   <- c(1.5, 2)

  xval <- xbart(x, z, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
                n.reps = n.reps,
                power = power,
                n.threads = 2L)
  
  expect_is(xval, "matrix")
  expect_equal(dim(xval), c(n.reps, length(power)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    power   = as.character(power)))
})


