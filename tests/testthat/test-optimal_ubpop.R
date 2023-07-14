context("test optimal_ubpop for correctness")
# test optimal_ubpop function
data(nydf, package = "smerc")
cases <- nydf$cases
pop <- nydf$population
coords <- nydf[, c("x", "y")]

estats <- suppressWarnings(optimal_ubpop(
  coords = coords,
  cases = cases,
  pop = pop,
  alpha = 1,
  nsim = 9,
  ubpop_seq = seq(0.05, 0.5, by = 0.05),
  min.cases = 0
))

# correct values of neg_lrt and gini_coef
# (created manually at some point)
neg_lrt <- c(
  -73.84, -70.72, -70.72, -70.72, -70.72, -70.72,
  -70.72, -70.72, -70.72, -70.72
)
gini_coef <- c(
  -0.23, -0.23, -0.23, -0.23, -0.23, -0.23, -0.23,
  -0.23, -0.23, -0.23
)

test_that("check accuracy for optimal_ubpop for NY data", {
  expect_equal(neg_lrt, round(estats$elbow_method$stats, 2))
  expect_equal(gini_coef, round(estats$gini_method$stats, 2))
  # check for correct upper bound
  expect_equal(0.1, estats$elbow_ubpop)
  expect_equal(0.1, estats$gini_ubpop)
})
