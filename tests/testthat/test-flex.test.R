set.seed(1)
coords = runif(8)

test_that("sanity checks for flex.test arguments", {
  expect_that(flex.test(coords), throws_error())
  coords = data.frame(x = runif(4), y = runif(4))
  cases = 1:3
  expect_that(flex.test(coords, cases = cases), throws_error())
  cases = 1:4
  expect_that(flex.test(coords, cases = as.factor(cases)), throws_error())
  pop = list(1:3)
  expect_that(flex.test(coords, cases = cases, pop = pop), throws_error())
  pop = list(1:4)
  expect_that(flex.test(coords, cases = cases, pop = pop), throws_error())
  pop = 1:4
  expect_that(flex.test(coords, cases = cases, pop = factor(pop)), throws_error())
  ex = 1:3
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex), throws_error())
  ex = 1:4
  alpha = -1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 1.1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = c(0.1, 0.3)
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 0.1
  nsim = 0
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim), 
              throws_error())
  nsim = c(10, 20)
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim), 
              throws_error())
  nsim = 10
  k = 0.5
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k), throws_error())
  k = c(1, 2)
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim,
                          k = k), throws_error())
  k = 2
  longlat = 1:2
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, longlat = longlat), throws_error())
  longlat = 1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, longlat = longlat), throws_error())
  longlat = FALSE
  parallel = 1:2
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, longlat = longlat, parallel = parallel), throws_error())
  parallel = 1
  expect_that(flex.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          k = k, longlat = longlat, parallel = parallel), throws_error())
  parallel = TRUE
  w = 1:4
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = matrix(1:4)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = diag(3)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = matrix(factor(diag(4)), nrow = 4)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w), throws_error())
  w = diag(4)
  expect_that(flex.test(coords, cases = cases, pop = pop, w = w, k = 10), throws_error())
})

set.seed(16)
data(nydf)
data(nyw)
coords = cbind(nydf$longitude, nydf$latitude)
cases = floor(nydf$cases)
pop = nydf$population
outp = flex.test(coords = coords, cases = cases, pop = pop, 
                  w = nyw, k = 10, nsim = 99, alpha = 1)

outb = flex.test(coords = coords, cases = cases, pop = pop, 
                 w = nyw, k = 10, nsim = 99, alpha = 1, 
                 type = "binomial")

# results taken from flex_test_ny_poisson_10nn_cartesian
test_that("check accuracy for flex.test poisson with FlexScan original for NY data", {
  expect_equal(sort(outp$clusters[[1]]$locids), 
               c(85, 86, 88, 89, 90, 92, 93))
  expect_equal(round(outp$clusters[[1]]$max_dist, 6), 0.245662)
  expect_equal(outp$clusters[[1]]$cases, 39)
  expect_equal(round(outp$clusters[[1]]$exp, 4), 16.3981)
  expect_equal(round(outp$clusters[[1]]$smr, 5), 2.37832)
  expect_equal(round(outp$clusters[[1]]$loglik, 4), 11.6713)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp$clusters[[1]]$pvalue, 0.01)
  
  expect_equal(sort(outp$clusters[[2]]$locids), 
               c(1, 2, 13, 15, 47, 49, 51))
  expect_equal(round(outp$clusters[[2]]$max_dist, 7),  0.0507453)
  expect_equal(outp$clusters[[2]]$cases, 31)
  expect_equal(round(outp$clusters[[2]]$exp, 4), 13.4462)
  expect_equal(round(outp$clusters[[2]]$smr, 5), 2.30548)
  expect_equal(round(outp$clusters[[2]]$loglik, 5), 8.62939)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp$clusters[[2]]$pvalue, 0.05)
  
  expect_equal(sort(outp$clusters[[9]]$locids), 
               c(102, 103, 106))
  expect_equal(round(outp$clusters[[9]]$max_dist, 6), 0.194769)
  expect_equal(outp$clusters[[9]]$cases, 11)
  expect_equal(round(outp$clusters[[9]]$exp, 5), 4.76234)
  expect_equal(round(outp$clusters[[9]]$smr, 5), 2.30979)
  expect_equal(round(outp$clusters[[9]]$loglik, 5), 3.00674)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp$clusters[[9]]$pvalue, 1)
})

# results taken from flex_test_ny_binomial_10nn_cartesian
test_that("check accuracy for flex.test binomial with FlexScan original for NY data", {
  expect_equal(sort(outb$clusters[[1]]$locids), 
               c(85, 86, 88, 89, 90, 92, 93))
  expect_equal(round(outb$clusters[[1]]$max_dist, 6), 0.245662)
  expect_equal(outb$clusters[[1]]$cases, 39)
  expect_equal(outb$clusters[[1]]$pop, 31420)
  expect_equal(round(outb$clusters[[1]]$smr, 5), 2.37832)
  expect_equal(round(outb$clusters[[1]]$loglik, 4), 11.6797)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[1]]$pvalue, 0.02)
  
  expect_equal(sort(outb$clusters[[2]]$locids), 
               c(1, 2, 13, 15, 47, 49, 51))
  expect_equal(round(outb$clusters[[2]]$max_dist, 7),  0.0507453)
  expect_equal(outb$clusters[[2]]$cases, 31)
  expect_equal(outb$clusters[[2]]$pop, 25764)
  expect_equal(round(outb$clusters[[2]]$smr, 5), 2.30548)
  expect_equal(round(outb$clusters[[2]]$loglik, 5), 8.63552)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[2]]$pvalue, 0.09)
  
  expect_equal(sort(outb$clusters[[9]]$locids), 
               c(102, 103, 106))
  expect_equal(round(outb$clusters[[9]]$max_dist, 6), 0.194769)
  expect_equal(outb$clusters[[9]]$cases, 11)
  expect_equal(outb$clusters[[9]]$pop, 9125)
  expect_equal(round(outb$clusters[[9]]$smr, 5), 2.30979)
  expect_equal(round(outb$clusters[[9]]$loglik, 5), 3.00889)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[9]]$pvalue, 1)
})

