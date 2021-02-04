fpath = system.file("testdata",  package = "smerc")
fname = paste(fpath, "/pruned_seq.rda", sep = "")
load(fname)

data(nydf, package = "smerc")
cases = nydf$cases
pop = nydf$population
coords = nydf[,c("x", "y")]
ex = sum(cases) / sum(pop) * pop
nsim = 499
alpha = 0.1
ubpop_seq = seq(0.01, 0.5, len = 50)
longlat = FALSE
cl = NULL
type = "poisson"
min.cases = 2
simdist = "multinomial"

set.seed(28) # necessary seed
pruned_seq_test = seq_scan_test(coords = coords,
                           cases = cases, pop = pop,
                           ex = ex,
                           nsim = nsim, alpha = alpha,
                           ubpop_seq = ubpop_seq,
                           longlat = longlat, cl = cl,
                           type = type,
                           min.cases = min.cases,
                           simdist = simdist)


context("check seq_scan_test with reference")
test_that("seq_scan_test matches ref", {
    expect_equal(pruned_seq_test, pruned_seq)
})
