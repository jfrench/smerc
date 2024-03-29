# 0.4.5
- Fix small bug in computing p-value for bn.test (unmodified).
- Add test for bn.test (modified).
# 1.0
- Add lots of methods (dc.test, fast.test, mlink.test, rflex.test, etc.)
- Dramatic restructure of many functions.
# 1.1
- Add cepp.test (and related).
- Fix some documentation.
# 1.2
- Update *.test functions to return 'smerc_cluster' class instead of 'scan' class.
- Create print, plot, and summary functions for smerc_cluster class.
Add ex argument to bn.test.
- Improve speed of flex.zones function.
# 1.3
- Create flex_zones and flex_test to use C++ code for speed.
- Add vignette regarding how to use package.
# 1.4
- Add optimal_ubpop for choosing a population upper bound.
# 1.5
- Add function to easily extract clusters.
# 1.6
- Add constant-risk version of Moran's I.
- Correct parallel processing bug (thanks to Insang Song).
# 1.7
- Improve efficiency of scan.test and elliptic.test by removing redundancies.
- Replace NEWS file with NEWS.md file
- Fix bug in scan.test when the population size of each region is identical. 
- Add nysf and nysp data sets, along with better descriptions.
# 1.8
- Add neast.
- Add gedist function to remove dependency on sp package.
- Add primes100k data set to sysdata.rda to remove dependency on randtoolbox package.
- Add nclusters function.
- Remove nysf and nypoly rda files from data. Replace with 
functions that create them from nysf. This is related to the point below.
- Remove sp from Depends and move to Suggests.
