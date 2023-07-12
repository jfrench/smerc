context("check precog.test accuracy")
set.seed(112)
library(smerc)

data("neast")
data("neastw")

# setup up arguments
coords <-  neast[, c("x", "y"), drop = TRUE]
pop <- neast$population
cases <- neast$cases

# precog
out <- precog.test(coords = coords, cases = cases,
                   pop = pop, w = neastw, nsim = 999)

# correct clusters
cl1 <-  c(88, 96, 83, 78, 89, 91, 157, 127, 77, 141, 84,
          210, 182, 205, 198, 213, 199, 194, 172, 216)
cl2 <- c(161, 163, 196, 202, 179, 221, 201, 183, 171,
         102, 112)
cl3 <- c(24, 23, 230, 13, 227, 16)

test_that("check accuracy for precog.test with neast", {
  expect_equal(out$clusters[[1]]$locids, cl1)
  expect_equal(round(out$clusters[[1]]$r, 2), 22412.31)
  expect_equal(out$clusters[[1]]$pop, 5504246)
  expect_equal(out$clusters[[1]]$cases, 12646)
  expect_equal(round(out$clusters[[1]]$exp, 2), 10984.75)
  expect_equal(round(out$clusters[[1]]$smr, 2), 1.15)
  expect_equal(round(out$clusters[[1]]$rr, 2), 1.19)
  expect_equal(round(out$clusters[[1]]$loglik, 4), 148.8356)
  expect_equal(round(out$clusters[[1]]$test_statistic, 4), 148.8356)
  expect_equal(out$clusters[[1]]$pvalue, 0.001)

  expect_equal(out$clusters[[2]]$locids, cl2)
  expect_equal(round(out$clusters[[2]]$r, 2), 29469.36)
  expect_equal(out$clusters[[2]]$pop, 1583988)
  expect_equal(out$clusters[[2]]$cases, 3778)
  expect_equal(round(out$clusters[[2]]$exp, 2), 3161.14)
  expect_equal(round(out$clusters[[2]]$smr, 2), 1.20)
  expect_equal(round(out$clusters[[2]]$rr, 2), 1.21)
  expect_equal(round(out$clusters[[2]]$loglik, 5), 60.03691)
  expect_equal(round(out$clusters[[2]]$test_statistic, 5), 60.03691)
  expect_equal(out$clusters[[2]]$pvalue, 0.001)

  expect_equal(out$clusters[[3]]$locids, cl3)
  expect_equal(round(out$clusters[[3]]$r, 2), 7464.83)
  expect_equal(out$clusters[[3]]$pop, 987552)
  expect_equal(out$clusters[[3]]$cases, 2255)
  expect_equal(round(out$clusters[[3]]$exp, 2), 1970.84)
  expect_equal(round(out$clusters[[3]]$smr, 2), 1.14)
  expect_equal(round(out$clusters[[3]]$rr, 2), 1.15)
  expect_equal(round(out$clusters[[3]]$loglik, 5), 20.27554)
  expect_equal(round(out$clusters[[3]]$test_statistic, 5), 20.27554)
  expect_equal(out$clusters[[3]]$pvalue, 0.002)
})


