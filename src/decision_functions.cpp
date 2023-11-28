#include "CASSIA.h"

// [[Rcpp::export]]
Rcpp::List plant_decision(double C_roots,
                          double N_roots,
                          double C_fungal,
                          double optimal_root_funga_biomass_ratio,
                          double N_allo,
                          double max_C_allocation_CASSIA) {

    /*
     * Carbon Allocation as in Mycofon!
     *
     * N_roots should be the N availability in the orginal model!
     * Thus the condition is changed as well the condition was original
     * N_av > 0.01 kg - N kg soil-1
     * and
     * N_allo < 0.5 (N_uptake_root + N_allo)
     * also
     * allomax * photosynthesis
     */

    double allomax;
    if (N_roots < 0.01) { // TODO: consider this again
      allomax = 1 - (1 - pow(exp(-50*N_roots), 3));
    } else {
      allomax = 0.2;
    }

    double allo; // Allocates all if it can meet it's own demand
    if (N_allo < 0.5*N_roots) { // TODO: consider this again
      allo = N_allo / (N_roots + N_allo);
    } else {
      allo = 1;
    }

    double C_allo = std::max(std::min(allomax * C_roots, allo * (C_roots * optimal_root_funga_biomass_ratio) - C_fungal), 0.0);

    /*
     *  Carbon Allocation in Franklin 2014
     *
     *  C_f is defined within a set of equations, but basically the excess after growth and turnover
     *
     *  C_f = Photosynthesis - G_p / y_p - turnover
     */

    double C_f = std::max(max_C_allocation_CASSIA, 0.0); // Value from CASSIA

    return(Rcpp::List::create(Rcpp::_["Mycofon_demand"] = 1,
                              Rcpp::_["Mycofon_allocation"] = C_allo,
                              Rcpp::_["Franklin_demand"] = 1,
                              Rcpp::_["Franklin_allocation"] = C_f));
}

// [[Rcpp::export]]
Rcpp::List myco_decision(double C_fungal,
                         double N_fungal,
                         double C_roots,
                         double N_roots,
                         double NC_fungal_opt,
                         double growth_C,
                         double growth_N) {

  /*
   * Nitrogen allocation from Mycofon!
   *
   * N_allo = N_max (1 - NC_root / NC_rootopt)
   *
   * N_max should be uptake rather than biomass!
   */

  double N_allo = std::max(N_fungal*(1 - (N_roots / C_roots) / NC_fungal_opt), 0.0);

  /*
   * Nitrogen allocation from Franklin 2014
   *
   * N_p = Uptake(B) - y_f C_f N:C
   *
   * Where the Uptake(B) is the biomass which is devoted to uptake
   *    TODO: mantle / ERM divide could be used here!
   *    TODO: Needs to be uptake rather than biomass though!
   */

  double growth = myco_growth(C_fungal, N_fungal, growth_C, growth_N)[2]; // Nitrogen used output
  double N_p = std::max(N_fungal - growth, 0.0);

  return(Rcpp::List::create(Rcpp::_["Mycofon_demand"] = 1,
                            Rcpp::_["Mycofon_allocation"] = N_allo,
                            Rcpp::_["Franklin_demand"] = 1,
                            Rcpp::_["Franklin_allocation"] = N_p));
}
