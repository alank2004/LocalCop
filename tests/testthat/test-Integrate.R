#--- Copula tests ------------------------------------------------------

## library(LocalCop)
## library(TMB)
## library(testthat)

# Generate Values
x <- runif(1, 0, 20)
y <- runif(1, 1, 3)
theta <- runif(1, 1, 3)

data <- list(model = "integral_function_test")
parameters <- list(x = x, y = y, theta = theta)

# TMB version of the function
integrated_tmb <- TMB::MakeADFun(data = data, parameters = parameters, silent = TRUE)

# R versions of function and derivative
integrated_r <- function(x, y, theta) {
  return(0.5 * x * x * y * theta)
}

integrated_grad_r <- function(x, y, theta) {
  return(matrix(c(x * y * theta, 0.5 * x * x * theta, 0.5 * x * x * y), nrow = 1))
}

# test
test_that("TMB function evaluation is equivalent to analytic R solution", {
  val_r <- integrated_r(x, y, theta)
  val_tmb <- integrated_tmb$fn(c(x, y, theta))
  expect_equal(val_r, val_tmb)
})

test_that("TMB function gradient evaluation is equivalent to analytic R solution", {
  grad_r <- integrated_grad_r(x, y, theta)
  grad_tmb <- integrated_tmb$gr(c(x, y, theta))
  expect_equal(grad_r, grad_tmb)
})
