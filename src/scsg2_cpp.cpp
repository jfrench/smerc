#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double hash_cz(std::vector<bool> &cz, const NumericVector &lprimes) {
  double pid = 0;
  unsigned int cz_size = cz.size();
  for(unsigned int i = 0; i < cz_size; i++) {
    if (cz[i] == true) {
      pid += lprimes[i];
    }
  }
  return pid;
}

std::list<std::vector<bool>> link_cz_nb(std::vector<bool> &cz,
                                           const IntegerVector &nb,
                                           std::unordered_set<double> &lz_hash_table,
                                           const NumericVector &lprimes) {
  unsigned int nb_size = nb.length();
  std::list<std::vector<bool>> lz;
  for(unsigned int i = 0; i < nb_size; i++) {
    cz[nb[i]] = true;
    if (lz_hash_table.insert(hash_cz(cz, lprimes)).second) {
      // insertion successful, add cz to list
      lz.emplace_back(cz);
    }
    cz[nb[i]] = false;
  }
  return lz;
}

IntegerVector colsums_sub_iv(IntegerMatrix &cw, IntegerVector sub) {
  unsigned int cw_ncol = cw.ncol();
  unsigned int sub_length = sub.length();
  IntegerVector msum(cw_ncol);
  for (unsigned int i = 0; i < cw_ncol; i++) {
    for (unsigned int j = 0; j < sub_length; j++) {
      msum[i] += cw(sub[j], i);
    }
  }
  return msum;
}

IntegerVector colsums_sub(IntegerMatrix &cw, std::vector<bool> sub) {
  unsigned int cw_ncol = cw.ncol();
  unsigned int sub_size = sub.size();
  IntegerVector msum(cw_ncol);
  for (unsigned int i = 0; i < cw_ncol; i++) {
    for (unsigned int j = 0; j < sub_size; j++) {
      if (sub[j]) {
        msum[i] += cw(j, i);
      }
    }
  }
  return msum;
}

IntegerVector add_biv(std::vector<bool> &lv,
                      IntegerVector iv) {
  unsigned int lv_size = lv.size();
  for (unsigned int i = 0; i < lv_size; i++) {
    if (lv[i] == true) {
      iv[i] += 1;
    }
  }
  return iv;
}

LogicalVector lmb(LogicalVector a, std::vector<bool> &b) {
  unsigned int a_length = a.length();
  for (unsigned int i = 0; i < a_length; i++) {
    if (b[i]) {
      a[i] = false;
    }
  }
  return a;
}

std::list<vector<bool>> csg2_cpp(std::vector<bool> &cz,
                                 IntegerVector &cnn,
                                 IntegerMatrix &cw,
                                 IntegerVector &s,
                                 std::unordered_set<double> &lz_hash_table,
                                 NumericVector &lprimes) {
  LogicalVector ab = (add_biv(cz, colsums_sub(cw, cz)) >= 1);
  return link_cz_nb(cz, s[lmb(ab, cz)], lz_hash_table, lprimes);
}

std::list<std::vector<bool>> lcsg2_cpp(std::list<std::vector<bool>> &lz,
                                       IntegerVector &cnn,
                                       IntegerMatrix &cw,
                                       IntegerVector &s,
                                       NumericVector &lprimes) {
  // initialize needed variables
  std::list<std::vector<bool>> new_lz;
  std::unordered_set<double> lz_hash_table;
  auto lz_end = lz.end();

  for (auto lit = lz.begin(); lit != lz_end; lit++) {
    // *lit has cz (current zone)
    // expand cz to lz and splice to end of new_lz
    new_lz.splice(new_lz.end(), csg2_cpp(*lit, cnn, cw, s, lz_hash_table, lprimes));
  }
  return new_lz;
}

IntegerMatrix sub_cnn(IntegerMatrix &x, IntegerVector &cnn) {

  // Determine the number of observations
  unsigned int nnn = cnn.size();

  // Create an output matrix
  IntegerMatrix out = no_init(nnn, nnn);

  // Loop through each column and copy the data.
  for(unsigned int i = 0; i < nnn; i++) {
    for (unsigned int j = 0; j < nnn; j++) {
      out(i, j) = x(cnn[i] - 1, cnn[j] - 1);
    }
  }
  return out;
}

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::export]]
std::list<std::list<std::vector<bool>>> scsg2_cpp(List nn,
                                              IntegerMatrix w,
                                              IntegerVector idx,
                                              unsigned int nlevel,
                                              NumericVector lprimes,
                                              bool verbose=false) {
  // declare int
  unsigned int i, clevel, nidx = idx.length();

  // current nearest neighbors, all region ids
  IntegerVector cnn, s;
  IntegerMatrix cw;

  // current list of zones
  std::list<std::vector<bool>> clz;
  // list of zones for idx i
  std::list<std::list<std::vector<bool>>> zidx;
  // list of zones
  std::list<std::list<std::vector<bool>>> z;

  // display progress
  Progress p(nidx, verbose);

  // loop over desired indices
  for (i = 0; i < nidx; i++) {
    p.increment();
    clevel = 1;
    // get current nn
    cnn = nn[idx[i] - 1];
    // get size of cnn vector
    unsigned int cnn_size = cnn.size();

    // create current w
    cw = sub_cnn(w, cnn);
    // create sequence for cnn
    s = seq(0, cnn_size - 1);

    // create current zone
    std::vector<bool> cz(cnn_size);
    cz[0] = true;

    // clear than create new clz for current zone
    clz.clear();
    clz.push_back(cz);
    // add clz to end of zidx
    zidx.clear();
    zidx.emplace_back(clz);

    // while there are still levels to consider
    while (clevel < nlevel) {
      // increment counter
      clevel++;
      // determine new set of zones from
      // previous list of candidate zones
      // and add to back list of zones
      zidx.emplace_back(lcsg2_cpp(zidx.back(), cnn, cw, s, lprimes));
      //if there are no more connected neighbors for
      // the current vector of zones, end expansion
      if (zidx.back().empty()) {
        clevel = nlevel;
      }
    }
    // clear clz to prepare for storage
    clz.clear();
    // move zidx indices to stacked list of bool vectors
    auto zidx_end = zidx.end();
    for (auto zidx_it = zidx.begin(); zidx_it != zidx_end; zidx_it++) {
      clz.splice(clz.end(), *zidx_it);
    }
    // move complete list of zones for region i to z
    z.emplace_back(std::move(clz));
  }
  return z;
}
