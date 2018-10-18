set.seed(15)
data(nydf)
data(nyw)
out = rflex.test(coords = cbind(nydf$longitude, nydf$latitude), 
                 cases = floor(nydf$cases), 
                 pop = nydf$population, 
                 w = nyw, k = 50,
                 nsim = 99, alpha = 1, longlat = FALSE,
                 alpha1 = 0.2)

test_that("check accuracy for rflex.test with FlexScan original for NY data", {
  expect_equal(sort(out$clusters[[1]]$locids), 
               c(1, 2, 13, 15, 27, 35, 37, 38, 43, 46, 47, 49, 51, 52, 53))
  expect_equal(round(out$clusters[[1]]$r, 5), 0.22483)
  expect_equal(out$clusters[[1]]$cases, 79)
  expect_equal(round(out$clusters[[1]]$exp, 4), 36.1025)
  expect_equal(round(out$clusters[[1]]$smr, 5), 2.18822)
  expect_equal(round(out$clusters[[1]]$loglik, 4), 20.8014)

  expect_equal(sort(out$clusters[[2]]$locids), 
               c(89, 90))
  expect_equal(round(out$clusters[[2]]$r, 6), 0.125172)
  expect_equal(out$clusters[[2]]$cases, 13)
  expect_equal(round(out$clusters[[2]]$exp, 5), 3.46124)
  expect_equal(round(out$clusters[[2]]$smr, 5), 3.75588)
  expect_equal(round(out$clusters[[2]]$loglik, 5), 7.74784)
  
  expect_equal(sort(out$clusters[[12]]$locids), 
               c(135, 146, 208, 210))
  expect_equal(round(out$clusters[[12]]$r, 6), 0.030934)
  expect_equal(out$clusters[[12]]$cases, 12)
  expect_equal(round(out$clusters[[12]]$exp, 5), 5.59582)
  expect_equal(round(out$clusters[[12]]$smr, 5), 2.14446)
  expect_equal(round(out$clusters[[12]]$loglik, 5), 2.78814)
})
