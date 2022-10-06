#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector stat_poisson_adj_cpp(NumericVector yin,
                                   NumericVector yout,
                                   NumericVector logein,
                                   NumericVector logeout,
                                   double a,
                                   NumericVector pen,
                                   unsigned int min_cases,
                                   bool return_max) {
  unsigned int yin_length = yin.length();
  double lrin, lrout;
  NumericVector tall(yin_length, 0);
  NumericVector tmax (1);

  // determine if there will be any problematic statistics
  for (unsigned int i = 0; i < yin_length; i++) {
    // compute statistic for good locations
    // yin > 0 and yin/ein > yout/ein
    if (yin[i] >= min_cases) {
      lrin = log(yin[i]) - logein[i];
      lrout = log(yout[i]) - logeout[i];
      if (lrin > lrout) {
        tall[i] = yin[i] * lrin + yout[i] * lrout;
      }
    }
  }

  // update test statistic if a > 0 for all shape values > 1
  if (a > 0) {
    // elliptical regions
    for (unsigned int j = 0; j < yin_length; j++) {
      if (pen[j] < 1) {
        tall[j] = tall[j] * pen[j];
      }
    }
  }
  if (return_max) {
    double dmax = max(tall);
    tmax[0] = dmax;
    return tmax;
  } else {
    return tall;
  }
}
