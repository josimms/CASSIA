#include "mycomodel.h"

// [[Rcpp::export]]
Rcpp::List myco_growth(double C_fungal, // per biomass
                       double N_fungal, // per biomass
                       double a,
                       double b)
{
  
  double carbon_used;
  double nitrogen_used;
  double out_1 = std::min(a*C_fungal, b*N_fungal); // The growth is just a experiential function, limited by the element
  double out_2 = std::max(out_1, 0.0); // No negative growth
  if (out_2 == a*C_fungal) {
    carbon_used = a*C_fungal;
    nitrogen_used = b*N_fungal*(a*C_fungal)/(b*N_fungal);
  } else if (out_2 == b*N_fungal) {
    carbon_used = a*C_fungal*(b*N_fungal)/(a*C_fungal);
    nitrogen_used = b*N_fungal;
  }
  
  return Rcpp::List::create(Rcpp::_["growth"] = out_2,
                            Rcpp::_["carbon_used"] = carbon_used,
                            Rcpp::_["nitrogen_used"] = nitrogen_used);
}