#--- Copula tests ------------------------------------------------------

## library(LocalCop)
## library(TMB)
## library(testthat)


# TMB version of the function
integrated_exp_tmb <- TMB::MakeADFun(data = list(model = "exponential_integral", lambda = 1, lower = 0, upper = 1), parameters = list(), silent = TRUE)

# R versions of function and derivative
integrated_exp_r <- function(lambda, lower, upper) {
  return(exp(-lambda * lower) - exp(-lambda * upper))
}

integrated_exp_grad_r <- function(lambda, lower, upper) {
  return(c(lambda * exp(-lambda * upper), -lambda * exp(-lambda * lower), upper * exp(-lambda * upper) - lower * exp(-lambda * lower)))
}

# test
lower <- runif(1)
upper <- lower + runif(1)
rate <- runif(1)

test_that("TMB function evaluation is equivalent to analytic R solution", {
  val_r <- integrated_exp_r(rate, lower, upper)
  val_tmb <- integrated_exp_tmb$fn(c(rate, lower, upper))
  expect_equal(val_r, val_tmb)
})

test_that("TMB function gradient evaluation is equivalent to analytic R solution", {
  grad_r <- integrated_exp_grad_r(rate, lower, upper)
  grad_tmb <- integrated_exp_tmb$gr(c(rate, lower, upper))
  expect_equal(grad_r, grad_tmb)
})
