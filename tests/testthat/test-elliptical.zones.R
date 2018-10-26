# preliminaries
set.seed(10)
data(nydf)
coords = nydf[,c("longitude", "latitude")]
pop = nydf$population
cases = nydf$cases
zones = smerc::elliptic.zones(coords, pop)

ty = sum(nydf$cases)
e = ty/sum(pop)*pop
ein = sapply(zones$zones, function(x) sum(e[x]))
eout = ty - ein
yin = sapply(zones$zones, function(x) sum(cases[x]))

tobs = scan.stat(yin, ein, eout, ty)

out = elliptic.test(coords, cases, pop, nsim = 0, a = 0)

context("check elliptic.zones")

test_that("compare results elliptic.zones and elliptic.test", {
  expect_equal(out$clusters[[1]]$test_statistic, max(tobs))
})

