#include "CASSIA.h"

// [[Rcpp::export]]
Rcpp::List myco_growth(double C_fungal,
                       double N_fungal,
                       double C_fungal_biomass,
                       double C_ecto,
                       double a,
                       double b,
                       double CN_ratio)
{

  double C_limitation = std::min(std::max(C_fungal / C_fungal_biomass, 0.0), 1.0); // TODO: percentage of storage
  double N_limitation = std::min(std::max(N_fungal / (C_fungal_biomass*CN_ratio), 0.0), 1.0); // TODO: percentage of storage

  double out_1 = C_ecto * C_limitation * N_limitation; // Limited by carbon
  double out_2 = std::max(out_1, 0.0); // No negative growth

  double carbon_used = out_2;
  double nitrogen_used = out_2 * CN_ratio;

  return Rcpp::List::create(Rcpp::_["growth"] = out_2,
                            Rcpp::_["carbon_used"] = carbon_used,
                            Rcpp::_["nitrogen_used"] = nitrogen_used);
}
