context("test optimal_ubpop for correctness")
# test optimal_ubpop function
data(nydf, package = "smerc")
cases = nydf$cases
pop = nydf$population
coords = nydf[,c("x", "y")]

set.seed(28)
estats = optimal_ubpop(coords = coords,
                       cases  = cases,
                       pop = pop,
                       alpha = 0.05,
                       nsim = 999,
                       ubpop_seq = seq(0.01, 0.5, len = 50),
                       min.cases = 0)
set.seed(28)
estats2 = optimal_ubpop2(coords = coords,
                       cases  = cases,
                       pop = pop,
                       alpha = 0.05,
                       nsim = 999,
                       ubpop_seq = seq(0.01, 0.5, len = 50),
                       min.cases = 0)

# all.equal(estats$neg_lrt, estats2$elbow_method$sig)
#
# all.equal(estats$sig_yin_seq, estats2$gini_method$sig_yin_seq)
# all.equal(estats$sig_ein_seq, estats2$gini_method$sig_ein_seq)
# all.equal(estats$gini_seq, estats2$gini_method$gini_seq)
# all.equal(estats$gini_seq * estats$all_sig, estats2$gini_method$gini_seq_sig)

# correct values of neg_lrt and gini_coef
# (created manually at some point)
neg_lrt =  c(-6.661,-7.115,-15.019,-15.019,-15.531,
             -15.779,-16.627,-17.252,-18.311,-21.03,
             -21.03,-21.03,-21.03,-21.03,-21.03,-21.03,
             -21.03,-21.03,-21.03,-21.03,-21.03,-21.03,
             -21.03,-21.03,-21.03,-21.03,-21.03,-21.03,
             -21.03,-21.03,-21.03,-21.03,-21.03,-21.03,
             -21.03,-21.03,-21.03,-21.03,-21.03,-21.03,
             -21.03,-21.03,-21.03,-21.03,-21.03,-21.03,
             -21.03,-21.03,-21.03,-21.03)
gini_coef = c(0,0,0.061,0.061,0.069,0.081,0.086,0.09,
              0.097,0.106,0.106,0.106,0.106,0.106,0.106,
              0.106,0.106,0.106,0.106,0.106,0.106,0.106,
              0.106,0.106,0.106,0.106,0.106,0.106,0.106,
              0.106,0.106,0.106,0.106,0.106,0.106,0.106,
              0.106,0.106,0.106,0.106,0.106,0.106,0.106,
              0.106,0.106,0.106,0.106,0.106,0.106,0.106)

test_that("check accuracy for optimal_ubpop for NY data", {
  expect_equal(neg_lrt, round(estats$neg_lrt, 3))
  expect_equal(neg_lrt, round(estats2$elbow_method$stats, 3))
  expect_equal(gini_coef, round(estats$gini_coef, 3))
  expect_equal(gini_coef, round(estats2$gini_method$stats, 3))
  # check for correct upper bound
  expect_equal(0.1, estats$elbow_ubpop)
  expect_equal(0.1, estats$gini_ubpop)
})
