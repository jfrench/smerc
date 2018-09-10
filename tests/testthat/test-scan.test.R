set.seed(1)
coords = runif(8)

context("check scan.test")

test_that("sanity checks for scan.test arguments", {
  expect_that(scan.test(coords), throws_error())
  coords = data.frame(x = runif(4), y = runif(4))
  cases = 1:3
  expect_that(scan.test(coords, cases = cases), throws_error())
  cases = 1:4
  expect_that(scan.test(coords, cases = as.factor(cases)), throws_error())
  pop = list(1:3)
  expect_that(scan.test(coords, cases = cases, pop = pop), throws_error())
  pop = list(1:4)
  expect_that(scan.test(coords, cases = cases, pop = pop), throws_error())
  pop = 1:4
  expect_that(scan.test(coords, cases = cases, pop = factor(pop)), throws_error())
  ex = 1:3
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex), throws_error())
  ex = 1:4
  alpha = -1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 1.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = c(0.1, 0.3)
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha), throws_error())
  alpha = 0.1
  nsim = c(10, 20)
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim), 
              throws_error())
  nsim = 10
  ubpop = -0.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop), throws_error())
  ubpop = 1.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim,
                          ubpop = ubpop), throws_error())
  ubpop = 1.1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop), throws_error())
  ubpop = c(0.1, 0.2)
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop), throws_error())
  ubpop = 0.5
  longlat = 1:2
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, longlat = longlat), throws_error())
  longlat = 1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, longlat = longlat), throws_error())
  longlat = FALSE
  parallel = 1:2
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, longlat = longlat, parallel = parallel), throws_error())
  parallel = 1
  expect_that(scan.test(coords, cases = cases, pop = pop, ex = ex, alpha = alpha, nsim = nsim, 
                          ubpop = ubpop, longlat = longlat, parallel = parallel), throws_error())
})

data(nydf)
out = scan.test(coords = cbind(nydf$longitude, nydf$latitude), 
                cases = floor(nydf$cases), pop = nydf$population, 
                longlat = TRUE, nsim = 49, alpha = .50)
# clusters from satscan.  it's not clear how satscan does 
# lon/lat distance.  they seem to match up very well
# with fields::rdist.earth, but sp::spDists should be more 
# accurate
cl1 = c(52, 50, 53, 38, 49, 48, 15, 39, 37, 1, 16, 44, 47, 
        40, 14, 2, 51, 13,
        43, 45, 17, 55, 11, 3, 12, 46, 36, 35, 54, 10, 5)
cl2 = c(88, 87, 92, 86, 89, 91, 93, 85, 90)
cl3 = c(113, 117, 116, 112, 220, 118, 115, 123, 
        124, 111, 114, 125, 219, 122, 126, 119)
diff1 = setdiff(out$clusters[[1]]$locids, cl1)
diff2 = setdiff(out$clusters[[2]]$locids, cl2)
diff3 = setdiff(out$clusters[[3]]$locids, cl3)

test_that("check accuracy for scan.test with SatScan for NY data", {
  expect_equal(length(diff1), 0)
  expect_equal(length(diff2), 0)
  expect_equal(length(diff3), 0)
  expect_equal(out$clusters[[1]]$pop, 119050)
  expect_equal(out$clusters[[1]]$cases, 106)
  expect_equal(round(out$clusters[[1]]$exp, 2), 62.13)
  expect_equal(round(out$clusters[[1]]$smr, 2), 1.71)
  expect_equal(round(out$clusters[[1]]$rr, 2), 1.87)
  expect_equal(round(out$clusters[[1]]$loglik, 2), 14.78)
  
  expect_that(out$clusters[[2]]$pop, equals(40696))
  expect_that(out$clusters[[2]]$cases, equals(42))
  expect_that(round(out$clusters[[2]]$exp, 2), equals(21.24))
  expect_that(round(out$clusters[[2]]$smr, 2), equals(1.98))
  expect_that(round(out$clusters[[2]]$rr, 2), equals(2.06))
  expect_that(round(out$clusters[[2]]$loglik, 2), equals(8.29))
  
  expect_that(out$clusters[[3]]$pop, equals(45667))
  expect_that(out$clusters[[3]]$cases, equals(44))
  expect_that(round(out$clusters[[3]]$exp, 2), equals(23.83))
  expect_that(round(out$clusters[[3]]$smr, 2), equals(1.85))
  expect_that(round(out$clusters[[3]]$rr, 2), equals(1.92))
  expect_that(round(out$clusters[[3]]$loglik, 2), equals(7.20))
})

out2 = scan.test(coords = cbind(nydf$longitude, nydf$latitude), 
                 cases = floor(nydf$cases), 
                 pop = nydf$population, 
                 longlat = FALSE, nsim = 0, alpha = 1)

test_that("check accuracy for scan.test with SatScan for NY data, cartesian", {
  locids1 = c(52, 50, 37, 49, 38, 48, 53, 39, 15, 47, 51, 1, 40, 44, 16, 2, 55, 36,
              14, 43, 13, 45, 3, 35, 17, 11, 12, 46)
  expect_equal(out2$clusters[[1]]$locids, locids1)
  expect_equal(round(out2$clusters[[1]]$r, 3), 0.087)
  expect_equal(out2$clusters[[1]]$pop, 111674)
  expect_equal(out2$clusters[[1]]$cases, 100)
  expect_equal(round(out2$clusters[[1]]$exp, 2), 58.28)
  expect_equal(round(out2$clusters[[1]]$smr, 2), 1.72)
  expect_equal(round(out2$clusters[[1]]$rr, 2), 1.87)
  expect_equal(round(out2$clusters[[1]]$loglik, 6), 14.083511)
  
  locids6 = c(166, 159, 167)
  expect_equal(out2$clusters[[6]]$locids, locids6)
  expect_equal(round(out2$clusters[[6]]$r, 3), 0.012)
  expect_equal(out2$clusters[[6]]$pop, 8839)
  expect_equal(out2$clusters[[6]]$cases, 11)
  expect_equal(round(out2$clusters[[6]]$exp, 2), 4.61)
  expect_equal(round(out2$clusters[[6]]$smr, 2), 2.38)
  expect_equal(round(out2$clusters[[6]]$rr, 2), 2.41)
  expect_equal(round(out2$clusters[[6]]$loglik, 6), 3.209485)
})


