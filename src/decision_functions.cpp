#include "CASSIA.h"

// [[Rcpp::export]]
Rcpp::List plant_decision(double C_roots_NonStruct,
                          double N_roots_NonStruct,
                          double C_fungal_NonStruct,
                          double optimal_root_funga_biomass_ratio) {

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

    /*
     * TODO: I am considering the non structural rather than actual values of C and N,
     * so the equations should be changed accordingly!
     */

    double allomax;
    if (N_roots_NonStruct < 0.01) { // TODO: consider this again
      allomax = 1 - (1 - pow(exp(-50*N_roots_NonStruct), 3));
    } else {
      allomax = 0.2;
    }

    double allo; // Allocates all if it can meet it's own demand
    if (N_roots_NonStruct < 0.5) { // TODO: consider this again
      allo = N_roots_NonStruct / (N_roots_NonStruct + 0.5);
    } else {
      allo = 1;
    }

    double temp = std::min(allomax * C_roots_NonStruct, allo * (C_roots_NonStruct * optimal_root_funga_biomass_ratio) - C_fungal_NonStruct);
    double C_allo = std::max(temp, 0.0);

    /*
     *  Carbon Allocation in Franklin 2014
     *
     *  C_f is defined within a set of equations, but basically the excess after growth and turnover
     *
     *  C_f = Photosynthesis - G_p / y_p - turnover
     */

    double C_f = std::max(C_roots_NonStruct, 0.0); // Value from CASSIA

    return(Rcpp::List::create(Rcpp::_["Mycofon_demand"] = 1,
                              Rcpp::_["Mycofon_allocation"] = C_allo,
                              Rcpp::_["Franklin_demand"] = 1,
                              Rcpp::_["Franklin_allocation"] = C_f));
}




// [[Rcpp::export]]
Rcpp::List myco_decision(double N_fungal_NonStruct,
                         double C_roots_NonStruct,
                         double N_roots_NonStruct,
                         double NC_fungal_opt) {

  /*
   * Nitrogen allocation from Mycofon!
   *
   * N_allo = N_max (1 - NC_root / NC_rootopt)
   *
   * N_max should be uptake rather than biomass!
   */

  double N_allo = std::max(N_fungal_NonStruct*(1 - (N_roots_NonStruct / C_roots_NonStruct) / NC_fungal_opt), 0.0);

  /*
   * Nitrogen allocation from Franklin 2014
   *
   * N_p = Uptake(B) - y_f C_f N:C
   *
   * Where the Uptake(B) is the biomass which is devoted to uptake
   *    TODO: mantle / ERM divide could be used here!
   *    TODO: Needs to be uptake rather than biomass though!
   */

  double N_p = std::max(N_fungal_NonStruct, 0.0);

  return(Rcpp::List::create(Rcpp::_["Mycofon_demand"] = 1,
                            Rcpp::_["Mycofon_allocation"] = N_allo,
                            Rcpp::_["Franklin_demand"] = 1,
                            Rcpp::_["Franklin_allocation"] = N_p));
}
