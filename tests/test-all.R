library(testthat)
library(smerc)
# only test in when long.double available
# don't have access to non-long.double machine for testing
if (capabilities("long.double")) {
  test_check("smerc")
}
