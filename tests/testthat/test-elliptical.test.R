set.seed(8)
data(nydf)
coords = nydf[,c("longitude", "latitude")]
pop = nydf$population
cases = floor(nydf$cases)
shape = c(1, 1.5, 2, 3, 4, 5)
nangle = c(1, 4, 6, 9, 12, 15)
ex = sum(cases)/sum(pop)*pop
ubpop = 0.5
cl = NULL
min.cases = 2

nsim = 0

out0 = elliptic.test(coords, cases, pop, nsim = nsim,
                    alpha = 0.5, a = 0)
out0.5 = elliptic.test(coords, cases, pop, nsim = nsim,
                     alpha = 0.5, a = 0.5)
out1 = elliptic.test(coords, cases, pop, nsim = nsim,
                       alpha = 1, a = 1)

locids0 = c(10, 5, 4, 18, 12, 11, 9, 6, 17, 3, 35, 7, 26, 13, 14, 2, 34, 16, 8,
           1, 15, 20, 47, 33, 27, 49, 48, 50, 24, 51, 25, 36, 19, 32, 52, 37,
           22, 23, 21, 39, 40, 38, 90, 28, 53, 55, 44, 83, 43, 89, 42, 45, 84,
           46, 54, 91, 41, 31, 250, 88, 87, 251, 237, 85, 238, 29, 240, 252,
           241, 172, 171, 249, 169, 227, 170, 168, 86, 92, 165, 163, 167, 166,
           164, 162, 159, 242, 93, 161, 225, 160, 158, 153, 154, 152, 151, 150,
           157, 148, 149, 236, 140, 142, 141, 155, 143, 139, 138, 147, 134, 156,
           133, 144, 125, 146, 132, 126, 137, 226, 121, 124, 122, 131, 127, 123,
           145, 136, 117, 120, 114, 128, 116, 118, 115, 135, 210, 239, 113, 111,
           129, 130, 119, 110, 112, 220, 219, 208, 224)
context("test elliptical.test")
test_that("check accuracy for elliptical.test a = 0", {
  expect_equal(locids0, 
               out0$clusters[[1]]$locids)
  expect_equal(0.26, 
               round(out0$clusters[[1]]$semiminor_axis, 2))
  expect_equal(1.02, 
               round(out0$clusters[[1]]$semimajor_axis, 2))
  expect_equal(-15, 90 - out0$clusters[[1]]$angle)
  expect_equal(4, out0$clusters[[1]]$shape)
  expect_equal(496712, out0$clusters[[1]]$pop)
  expect_equal(337, out0$clusters[[1]]$cases)
  expect_equal(259.23, round(out0$clusters[[1]]$ex, 2))
  expect_equal(1.3, round(out0$clusters[[1]]$smr, 2))
  expect_equal(1.77, round(out0$clusters[[1]]$rr, 2))
  expect_equal(22.034279, round(out0$clusters[[1]]$test_statistic, 6))
  expect_equal(22.034279, round(out0$clusters[[1]]$loglikrat, 6))
})

test_that("check accuracy for elliptical.test a = 0.5", {
  expect_equal(locids0, 
               out0.5$clusters[[1]]$locids)
  expect_equal(0.26, 
               round(out0.5$clusters[[1]]$semiminor_axis, 2))
  expect_equal(1.02, 
               round(out0.5$clusters[[1]]$semimajor_axis, 2))
  expect_equal(-15, 90 - out0.5$clusters[[1]]$angle)
  expect_equal(4, out0.5$clusters[[1]]$shape)
  expect_equal(496712, out0.5$clusters[[1]]$pop)
  expect_equal(337, out0.5$clusters[[1]]$cases)
  expect_equal(259.23, round(out0.5$clusters[[1]]$ex, 2))
  expect_equal(1.3, round(out0.5$clusters[[1]]$smr, 2))
  expect_equal(1.77, round(out0.5$clusters[[1]]$rr, 2))
  expect_equal(17.627423, round(out0.5$clusters[[1]]$test_statistic, 6))
  expect_equal(22.034279, round(out0.5$clusters[[1]]$loglikrat, 6))
})

locids1.1 = c(52, 53, 50, 38, 15, 49, 16, 44, 1, 48, 14, 39, 13, 2, 45, 17, 37, 43,
              11, 47, 46, 12, 40, 54, 3, 10, 51, 5, 18, 9, 41, 35)
locids1.6 = c(266, 281, 265)
locids1.7 = c(106, 103, 102, 77, 230)

test_that("check accuracy for elliptical.test a = 1", {
  expect_equal(locids1.1, 
               out1$clusters[[1]]$locids)
  expect_equal(0.056, 
               round(out1$clusters[[1]]$semiminor_axis, 3))
  expect_equal(0.11, 
               round(out1$clusters[[1]]$semimajor_axis, 2))
  expect_equal(90 + 90, out1$clusters[[1]]$angle)
  expect_equal(2, out1$clusters[[1]]$shape)
  expect_equal(124896, out1$clusters[[1]]$pop)
  expect_equal(114, out1$clusters[[1]]$cases)
  expect_equal(65.18, round(out1$clusters[[1]]$ex, 2))
  expect_equal(1.75, round(out1$clusters[[1]]$smr, 2))
  expect_equal(1.94, round(out1$clusters[[1]]$rr, 2))
  expect_equal(15.504490, round(out1$clusters[[1]]$test_statistic, 6))
  expect_equal(17.442551, round(out1$clusters[[1]]$loglikrat, 6))
  
  expect_equal(locids1.6, 
               out1$clusters[[6]]$locids)
  expect_equal(0.017, 
               round(out1$clusters[[6]]$semiminor_axis, 3))
  expect_equal(0.026, 
               round(out1$clusters[[6]]$semimajor_axis, 3))
  expect_equal(90 + 90, out1$clusters[[6]]$angle)
  expect_equal(1.5, out1$clusters[[6]]$shape)
  expect_equal(8063, out1$clusters[[6]]$pop)
  expect_equal(11, out1$clusters[[6]]$cases)
  expect_equal(4.21, round(out1$clusters[[6]]$ex, 2))
  expect_equal(2.61, round(out1$clusters[[6]]$smr, 2))
  expect_equal(2.65, round(out1$clusters[[6]]$rr, 2))
  expect_equal(3.667327, round(out1$clusters[[6]]$test_statistic, 6))
  expect_equal(3.820132, round(out1$clusters[[6]]$loglikrat, 6))
  
  expect_equal(locids1.7, 
               out1$clusters[[7]]$locids)
  expect_equal(0.087, 
               round(out1$clusters[[7]]$semiminor_axis, 3))
  expect_equal(0.35, 
               round(out1$clusters[[7]]$semimajor_axis, 2))
  expect_equal(90 + 45, out1$clusters[[7]]$angle)
  expect_equal(4.00, out1$clusters[[7]]$shape)
  expect_equal(15576, out1$clusters[[7]]$pop)
  expect_equal(19, out1$clusters[[7]]$cases)
  expect_equal(8.13, round(out1$clusters[[7]]$ex, 2))
  expect_equal(2.34, round(out1$clusters[[7]]$smr, 2))
  expect_equal(2.38, round(out1$clusters[[7]]$rr, 2))
  expect_equal(3.436309, round(out1$clusters[[7]]$test_statistic, 6))
  expect_equal(5.369233, round(out1$clusters[[7]]$loglikrat, 6))
})

