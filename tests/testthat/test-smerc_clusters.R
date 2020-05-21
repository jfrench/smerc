# test updated methods

data(nydf)
data(nyw)
coords = with(nydf, cbind(longitude, latitude))
set.seed(1)
out = uls.test(coords = coords, cases = floor(nydf$cases),
               pop = nydf$pop, w = nyw,
               alpha = 0.05, longlat = TRUE,
               nsim = 19, ubpop = 0.5)
set.seed(1)
out2 = uls.test2(coords = coords, cases = floor(nydf$cases),
               pop = nydf$pop, w = nyw,
               alpha = 0.05, longlat = TRUE,
               nsim = 19, ubpop = 0.5)

set.seed(1)
out3 = rflex.test(coords = coords, cases = floor(nydf$cases),
                  w = nyw, k = 15,
                  pop = nydf$pop, nsim = 19,
                  alpha = 0.1, longlat = TRUE)

set.seed(1)
out4 = rflex.test2(coords = coords, cases = floor(nydf$cases),
                   w = nyw, k = 15,
                   pop = nydf$pop, nsim = 19,
                   alpha = 0.1, longlat = TRUE)

set.seed(1)
out5 = mlink.test(coords = coords, cases = floor(nydf$cases),
                  pop = nydf$pop, w = nyw,
                  alpha = 0.12, longlat = TRUE,
                  nsim = 19, ubpop = 0.1, ubd = 0.2)

set.seed(1)
out6 = mlink.test2(coords = coords, cases = floor(nydf$cases),
                   pop = nydf$pop, w = nyw,
                   alpha = 0.12, longlat = TRUE,
                   nsim = 19, ubpop = 0.1, ubd = 0.2)

set.seed(1)
out7 = flex.test(coords = coords, cases = floor(nydf$cases),
                 w = nyw, k = 10,
                 pop = nydf$pop, nsim = 19,
                 alpha = 0.12, longlat = TRUE)
set.seed(1)
out8 = flex.test2(coords = coords, cases = floor(nydf$cases),
                  w = nyw, k = 10,
                  pop = nydf$pop, nsim = 19,
                  alpha = 0.12, longlat = TRUE)

set.seed(1)
out9 = fast.test(coords = coords, cases = floor(nydf$cases),
                 pop = nydf$pop,
                 alpha = 0.1, longlat = TRUE,
                 nsim = 19, ubpop = 0.5)
set.seed(1)
out10 = fast.test2(coords = coords, cases = floor(nydf$cases),
                   pop = nydf$pop,
                   alpha = 0.1, longlat = TRUE,
                   nsim = 19, ubpop = 0.5)

set.seed(1)
out11 = mlf.test(coords = coords, cases = floor(nydf$cases),
                 pop = nydf$pop, w = nyw,
                 alpha = 0.12, longlat = TRUE,
                 nsim = 19, ubpop = 0.1, ubd = 0.5)
set.seed(1)
out12 = mlf.test2(coords = coords, cases = floor(nydf$cases),
                  pop = nydf$pop, w = nyw,
                  alpha = 0.12, longlat = TRUE,
                  nsim = 19, ubpop = 0.1, ubd = 0.5)

set.seed(1)
out13 = elliptic.test(coords = coords,
                      cases = floor(nydf$cases),
                      pop = nydf$pop, ubpop = 0.1,
                      nsim = 19,
                      alpha = 0.12)
set.seed(1)
out14 = elliptic.test2(coords = coords,
                       cases = floor(nydf$cases),
                       pop = nydf$pop, ubpop = 0.1,
                       nsim = 19,
                       alpha = 0.12)

set.seed(1)
edmst_test_ref = edmst.test(coords = coords, cases = floor(nydf$cases),
                            pop = nydf$pop, w = nyw,
                            alpha = 0.12, longlat = TRUE,
                            nsim = 19, ubpop = 0.1, ubd = 0.2)
set.seed(1)
edmst_test_check = edmst.test2(coords = coords, cases = floor(nydf$cases),
                              pop = nydf$pop, w = nyw,
                              alpha = 0.12, longlat = TRUE,
                              nsim = 19, ubpop = 0.1, ubd = 0.2)

set.seed(1)
dmst_test_ref = dmst.test(coords = coords, cases = floor(nydf$cases),
                            pop = nydf$pop, w = nyw,
                            alpha = 0.12, longlat = TRUE,
                            nsim = 19, ubpop = 0.1, ubd = 0.2)
set.seed(1)
dmst_test_check = dmst.test2(coords = coords, cases = floor(nydf$cases),
                               pop = nydf$pop, w = nyw,
                               alpha = 0.12, longlat = TRUE,
                               nsim = 19, ubpop = 0.1, ubd = 0.2)

set.seed(1)
dc_test_ref = dc.test(coords = coords, cases = floor(nydf$cases),
                          pop = nydf$pop, w = nyw,
                          alpha = 0.12, longlat = TRUE,
                          nsim = 19, ubpop = 0.1, ubd = 0.2)
set.seed(1)
dc_test_check = dc.test2(coords = coords, cases = floor(nydf$cases),
                             pop = nydf$pop, w = nyw,
                             alpha = 0.12, longlat = TRUE,
                             nsim = 19, ubpop = 0.1, ubd = 0.2)

library(testthat)
test_that("uls.test and uls2 match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(out$clusters[[i]]$locids, out2$clusters[[i]]$locids)
    expect_equal(out$clusters[[i]]$coords, out2$clusters[[i]]$centroid)
    expect_equal(out$clusters[[i]]$r, out2$clusters[[i]]$r)
    expect_equal(out$clusters[[i]]$max_dist, out2$clusters[[i]]$max_dist)
    expect_equal(out$clusters[[i]]$pop, out2$clusters[[i]]$pop)
    expect_equal(out$clusters[[i]]$cases, out2$clusters[[i]]$cases)
    expect_equal(out$clusters[[i]]$ex, out2$clusters[[i]]$expected)
    expect_equal(out$clusters[[i]]$smr, out2$clusters[[i]]$smr)
    expect_equal(out$clusters[[i]]$rr, out2$clusters[[i]]$rr)
    expect_equal(out$clusters[[i]]$loglikrat, out2$clusters[[i]]$loglikrat)
    expect_equal(out$clusters[[i]]$test_statistic, out2$clusters[[i]]$test_statistic)
    expect_equal(out$clusters[[i]]$pvalue, out2$clusters[[i]]$pvalue)
    expect_equal(out$clusters[[i]]$w, out2$clusters[[i]]$w)
  }
})



test_that("rflex.test and rflex2 match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(out3$clusters[[i]]$locids, out4$clusters[[i]]$locids)
    expect_equal(out3$clusters[[i]]$coords, out4$clusters[[i]]$centroid)
    expect_equal(out3$clusters[[i]]$r, out4$clusters[[i]]$r)
    expect_equal(out3$clusters[[i]]$max_dist, out4$clusters[[i]]$max_dist)
    expect_equal(out3$clusters[[i]]$pop, out4$clusters[[i]]$pop)
    expect_equal(out3$clusters[[i]]$cases, out4$clusters[[i]]$cases)
    expect_equal(out3$clusters[[i]]$ex, out4$clusters[[i]]$expected)
    expect_equal(out3$clusters[[i]]$smr, out4$clusters[[i]]$smr)
    expect_equal(out3$clusters[[i]]$rr, out4$clusters[[i]]$rr)
    expect_equal(out3$clusters[[i]]$loglikrat, out4$clusters[[i]]$loglikrat)
    expect_equal(out3$clusters[[i]]$test_statistic, out4$clusters[[i]]$test_statistic)
    expect_equal(out3$clusters[[i]]$pvalue, out4$clusters[[i]]$pvalue)
    expect_equal(out3$clusters[[i]]$w, out4$clusters[[i]]$w)
  }
})


test_that("mlink.test and mlink2 match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(out5$clusters[[i]]$locids, out6$clusters[[i]]$locids)
    expect_equal(out5$clusters[[i]]$coords, out6$clusters[[i]]$centroid)
    expect_equal(out5$clusters[[i]]$r, out6$clusters[[i]]$r)
    expect_equal(out5$clusters[[i]]$max_dist, out6$clusters[[i]]$max_dist)
    expect_equal(out5$clusters[[i]]$pop, out6$clusters[[i]]$pop)
    expect_equal(out5$clusters[[i]]$cases, out6$clusters[[i]]$cases)
    expect_equal(out5$clusters[[i]]$ex, out6$clusters[[i]]$expected)
    expect_equal(out5$clusters[[i]]$smr, out6$clusters[[i]]$smr)
    expect_equal(out5$clusters[[i]]$rr, out6$clusters[[i]]$rr)
    expect_equal(out5$clusters[[i]]$loglikrat, out6$clusters[[i]]$loglikrat)
    expect_equal(out5$clusters[[i]]$test_statistic, out6$clusters[[i]]$test_statistic)
    expect_equal(out5$clusters[[i]]$pvalue, out6$clusters[[i]]$pvalue)
    expect_equal(out5$clusters[[i]]$w, out6$clusters[[i]]$w)
  }
})


test_that("flex.test and flex2 match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(out7$clusters[[i]]$locids, out8$clusters[[i]]$locids)
    expect_equal(out7$clusters[[i]]$coords, out8$clusters[[i]]$centroid)
    expect_equal(out7$clusters[[i]]$r, out8$clusters[[i]]$r)
    expect_equal(out7$clusters[[i]]$max_dist, out8$clusters[[i]]$max_dist)
    expect_equal(out7$clusters[[i]]$pop, out8$clusters[[i]]$pop)
    expect_equal(out7$clusters[[i]]$cases, out8$clusters[[i]]$cases)
    expect_equal(out7$clusters[[i]]$ex, out8$clusters[[i]]$expected)
    expect_equal(out7$clusters[[i]]$smr, out8$clusters[[i]]$smr)
    expect_equal(out7$clusters[[i]]$rr, out8$clusters[[i]]$rr)
    expect_equal(out7$clusters[[i]]$loglikrat, out8$clusters[[i]]$loglikrat)
    expect_equal(out7$clusters[[i]]$test_statistic, out8$clusters[[i]]$test_statistic)
    expect_equal(out7$clusters[[i]]$pvalue, out8$clusters[[i]]$pvalue)
    expect_equal(out7$clusters[[i]]$w, out8$clusters[[i]]$w)
  }
})


test_that("fast.test and fast2 match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(out9$clusters[[i]]$locids, out10$clusters[[i]]$locids)
    expect_equal(out9$clusters[[i]]$coords, out10$clusters[[i]]$centroid)
    expect_equal(out9$clusters[[i]]$r, out10$clusters[[i]]$r)
    expect_equal(out9$clusters[[i]]$max_dist, out10$clusters[[i]]$max_dist)
    expect_equal(out9$clusters[[i]]$pop, out10$clusters[[i]]$pop)
    expect_equal(out9$clusters[[i]]$cases, out10$clusters[[i]]$cases)
    expect_equal(out9$clusters[[i]]$ex, out10$clusters[[i]]$expected)
    expect_equal(out9$clusters[[i]]$smr, out10$clusters[[i]]$smr)
    expect_equal(out9$clusters[[i]]$rr, out10$clusters[[i]]$rr)
    expect_equal(out9$clusters[[i]]$loglikrat, out10$clusters[[i]]$loglikrat)
    expect_equal(out9$clusters[[i]]$test_statistic, out10$clusters[[i]]$test_statistic)
    expect_equal(out9$clusters[[i]]$pvalue, out10$clusters[[i]]$pvalue)
    expect_equal(out9$clusters[[i]]$w, out10$clusters[[i]]$w)
  }
})


test_that("mlf.test and mlf2 match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(out11$clusters[[i]]$locids, out12$clusters[[i]]$locids)
    # expect_equal(out11$clusters[[i]]$coords, out12$clusters[[i]]$centroid)
    # expect_equal(out11$clusters[[i]]$r, out12$clusters[[i]]$r)
    # expect_equal(out11$clusters[[i]]$max_dist, out12$clusters[[i]]$max_dist)
    expect_equal(out11$clusters[[i]]$pop, out12$clusters[[i]]$pop)
    expect_equal(out11$clusters[[i]]$cases, out12$clusters[[i]]$cases)
    expect_equal(out11$clusters[[i]]$ex, out12$clusters[[i]]$expected)
    expect_equal(out11$clusters[[i]]$smr, out12$clusters[[i]]$smr)
    expect_equal(out11$clusters[[i]]$rr, out12$clusters[[i]]$rr)
    expect_equal(out11$clusters[[i]]$loglikrat, out12$clusters[[i]]$loglikrat)
    # expect_equal(out11$clusters[[i]]$test_statistic, out12$clusters[[i]]$test_statistic)
    expect_equal(out11$clusters[[i]]$pvalue, out12$clusters[[i]]$pvalue)
    expect_equal(out11$clusters[[i]]$w, out12$clusters[[i]]$w)
  }
})


test_that("edmst_test_ref and edmst.test match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(edmst_test_ref$clusters[[i]]$locids, edmst_test_check$clusters[[i]]$locids)
    expect_equal(edmst_test_ref$clusters[[i]]$coords, edmst_test_check$clusters[[i]]$centroid)
    expect_equal(edmst_test_ref$clusters[[i]]$r, edmst_test_check$clusters[[i]]$r)
    expect_equal(edmst_test_ref$clusters[[i]]$max_dist, edmst_test_check$clusters[[i]]$max_dist)
    expect_equal(edmst_test_ref$clusters[[i]]$pop, edmst_test_check$clusters[[i]]$pop)
    expect_equal(edmst_test_ref$clusters[[i]]$cases, edmst_test_check$clusters[[i]]$cases)
    expect_equal(edmst_test_ref$clusters[[i]]$ex, edmst_test_check$clusters[[i]]$expected)
    expect_equal(edmst_test_ref$clusters[[i]]$smr, edmst_test_check$clusters[[i]]$smr)
    expect_equal(edmst_test_ref$clusters[[i]]$rr, edmst_test_check$clusters[[i]]$rr)
    expect_equal(edmst_test_ref$clusters[[i]]$loglikrat, edmst_test_check$clusters[[i]]$loglikrat)
    expect_equal(edmst_test_ref$clusters[[i]]$test_statistic, edmst_test_check$clusters[[i]]$test_statistic)
    expect_equal(edmst_test_ref$clusters[[i]]$pvalue, edmst_test_check$clusters[[i]]$pvalue)
    expect_equal(edmst_test_ref$clusters[[i]]$w, edmst_test_check$clusters[[i]]$w)
  }
})

test_that("dmst_test_ref and dmst.test match", {
  for(i in seq_along(out$clusters)) {
    expect_equal(dmst_test_ref$clusters[[i]]$locids, dmst_test_check$clusters[[i]]$locids)
    expect_equal(dmst_test_ref$clusters[[i]]$coords, dmst_test_check$clusters[[i]]$centroid)
    expect_equal(dmst_test_ref$clusters[[i]]$r, dmst_test_check$clusters[[i]]$r)
    expect_equal(dmst_test_ref$clusters[[i]]$max_dist, dmst_test_check$clusters[[i]]$max_dist)
    expect_equal(dmst_test_ref$clusters[[i]]$pop, dmst_test_check$clusters[[i]]$pop)
    expect_equal(dmst_test_ref$clusters[[i]]$cases, dmst_test_check$clusters[[i]]$cases)
    expect_equal(dmst_test_ref$clusters[[i]]$ex, dmst_test_check$clusters[[i]]$expected)
    expect_equal(dmst_test_ref$clusters[[i]]$smr, dmst_test_check$clusters[[i]]$smr)
    expect_equal(dmst_test_ref$clusters[[i]]$rr, dmst_test_check$clusters[[i]]$rr)
    expect_equal(dmst_test_ref$clusters[[i]]$loglikrat, dmst_test_check$clusters[[i]]$loglikrat)
    expect_equal(dmst_test_ref$clusters[[i]]$test_statistic, dmst_test_check$clusters[[i]]$test_statistic)
    expect_equal(dmst_test_ref$clusters[[i]]$pvalue, dmst_test_check$clusters[[i]]$pvalue)
    expect_equal(dmst_test_ref$clusters[[i]]$w, dmst_test_check$clusters[[i]]$w)
  }
})

test_that("dc_test_ref and dc.test match", {
  for (i in seq_along(out$clusters)) {
    expect_equal(dc_test_ref$clusters[[i]]$locids, dc_test_check$clusters[[i]]$locids)
    expect_equal(dc_test_ref$clusters[[i]]$coords, dc_test_check$clusters[[i]]$centroid)
    expect_equal(dc_test_ref$clusters[[i]]$r, dc_test_check$clusters[[i]]$r)
    expect_equal(dc_test_ref$clusters[[i]]$max_dist, dc_test_check$clusters[[i]]$max_dist)
    expect_equal(dc_test_ref$clusters[[i]]$pop, dc_test_check$clusters[[i]]$pop)
    expect_equal(dc_test_ref$clusters[[i]]$cases, dc_test_check$clusters[[i]]$cases)
    expect_equal(dc_test_ref$clusters[[i]]$ex, dc_test_check$clusters[[i]]$expected)
    expect_equal(dc_test_ref$clusters[[i]]$smr, dc_test_check$clusters[[i]]$smr)
    expect_equal(dc_test_ref$clusters[[i]]$rr, dc_test_check$clusters[[i]]$rr)
    expect_equal(dc_test_ref$clusters[[i]]$loglikrat, dc_test_check$clusters[[i]]$loglikrat)
    expect_equal(dc_test_ref$clusters[[i]]$test_statistic, dc_test_check$clusters[[i]]$test_statistic)
    expect_equal(dc_test_ref$clusters[[i]]$pvalue, dc_test_check$clusters[[i]]$pvalue)
    expect_equal(dc_test_ref$clusters[[i]]$w, dc_test_check$clusters[[i]]$w)
  }
})




