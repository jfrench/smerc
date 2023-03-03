context("check scan.sim accuracy for different distributions")
set.seed(2)
nsim <- 499
data(nydf)
coords <- nydf[, c("x", "y")]
nn <- nnpop(as.matrix(dist(coords)), pop = nydf$pop, ubpop = 0.1)
cases <- floor(nydf$cases)
ty <- sum(cases)
e <- ty / sum(nydf$population) * nydf$population
ein <- nn.cumsum(nn, e)
eout <- ty - ein
tpop <- sum(nydf$population)
popin <- nn.cumsum(nn, nydf$population)
popout <- tpop - popin
sa <- scan.sim.adj(nsim, nn,
  ty = ty, ex = e, type = "poisson",
  logein = log(ein), logeout = log(eout), simdist = "multinomial",
  pop = nydf$pop
)

sb <- scan.sim.adj(nsim, nn,
  ty = ty, ex = e, type = "poisson",
  logein = log(ein), logeout = log(eout),
  simdist = "poisson",
  pop = nydf$pop
)

sc <- scan.sim.adj(nsim, nn,
  ty = ty, ex = e, type = "binomial",
  logein = log(eout), logeout = log(eout), simdist = "binomial",
  tpop = tpop, popin = popin, logpopin = log(popin),
  popout = tpop - popin, logpopout = log(popout),
  pop = nydf$pop
)
summa <- summary(sa)
summb <- summary(sb)
summc <- summary(sc)

test_that("check accuracy for scan.sim", {
  expect_true(round(summa[2], 1) - round(summb[2], 1) <= 0.1)
  expect_true(round(summb[2], 1) - round(summc[2], 1) <= 0.1)
  expect_true(round(summa[3], 1) - round(summb[3], 1) <= 0.1)
  expect_true(round(summb[3], 1) - round(summc[3], 1) <= 0.1)
  expect_true(round(summa[4], 1) - round(summb[4], 1) <= 0.1)
  expect_true(round(summb[4], 1) - round(summc[4], 1) <= 0.1)
})
