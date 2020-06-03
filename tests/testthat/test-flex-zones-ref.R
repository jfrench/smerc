fpath = system.file("testdata",  package = "smerc")
fname = paste(fpath, "/flex10_zones_ref.rda", sep = "")
load(fname)

data(nydf)
coords = as.matrix(nydf[, c("x", "y")])
data(nyw)
flex10_zones = flex.zones(coords, nyw, k = 10)

context("check flex10_zones_ref w/ flex10_zones")
test_that("flex.zones result matches reference", {
  expect_equal(flex10_zones_ref, flex10_zones)
})
