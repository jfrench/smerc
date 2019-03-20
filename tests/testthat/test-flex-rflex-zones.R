# preliminaries
data(nydf)
data(nyw)
coords = nydf[, c("longitude", "latitude")]
pop = nydf$population
cases = nydf$cases

# construct zones
fzones = flex.zones(coords, w = nyw, k = 3, longlat = TRUE)
fzones2 = flex.zones(coords, w = nyw, k = 3, longlat = TRUE,
                     verbose = TRUE)

nn = knn(coords, longlat = TRUE, k = 10)
ex = pop * sum(cases) / sum(pop)
rzones = rflex.zones(nn, w = nyw, cases = floor(cases), ex = ex)
rzones2 = rflex.zones(nn, w = nyw, cases = floor(cases),
                      ex = ex, loop = TRUE)

context("check loop flex, rflex.zones")
test_that("compare flex.zones and rflex.zones w/ and w/o loop", {
  expect_equal(fzones, fzones2)
  expect_equal(rzones, rzones2)
})
