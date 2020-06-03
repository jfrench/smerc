# save output
fpath = system.file("testdata",  package = "smerc")

data(nydf)
coords = as.matrix(nydf[, c("x", "y")])
data(nyw)
flex10_zones_ref = flex.zones(coords, nyw, k = 10)

fname = paste(fpath, "/flex10_zones_ref.rda", sep = "")
save(flex10_zones_ref,
     compress = "bzip2",
     file = fname,
     version = 2)
