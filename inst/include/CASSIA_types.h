#include "function_structures.h"
// includes iostream, vector
#include <RcppCommon.h>

#ifndef CASSIA_TYPES_H
#define CASSIA_TYPES_H

namespace Rcpp {
template <>
SEXP wrap(const carbo_balance& x);
}

#endif

#include <Rcpp.h>

#ifndef CASSIA_TYPES_H
#define CASSIA_TYPES_H

namespace Rcpp {
template <>
SEXP wrap(const carbo_balance& x) {
  return Rcpp::List::create(Rcpp::Named("sugar_needles") = x.sugar.needles,
                            Rcpp::Named("sugar_phloem") = x.sugar.phloem,
                            Rcpp::Named("sugar_xylem_sh") = x.sugar.xylem_sh,
                            Rcpp::Named("sugar_xylem_st") = x.sugar.xylem_st,
                            Rcpp::Named("sugar_roots") = x.sugar.roots,
                            Rcpp::Named("starch_needles") = x.starch.needles,
                            Rcpp::Named("starch_phloem") = x.starch.phloem,
                            Rcpp::Named("starch_xylem_sh") = x.starch.xylem_sh,
                            Rcpp::Named("starch_xylem_st") = x.starch.xylem_st,
                            Rcpp::Named("starch_roots") = x.starch.roots);
}
}

#endif
