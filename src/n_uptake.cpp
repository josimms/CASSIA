#include "CASSIA.h"


// Input R values in the correct structure - N_balence
N_balence vector_to_N_balence(std::vector<double> input) {
  N_balence params;

  params.NH4 = input[0];
  params.NO3 = input[1];
  params.Norg = input[2];
  params.C = input[3];
  params.Norg_FOM = input[4];

  return(params);
}


// Input R values in the correct structure - N_balence
N_balence list_to_N_balence(Rcpp::List input) {
  N_balence params;

  params.NH4 = input[0];
  params.NO3 = input[1];
  params.Norg = input[2];
  params.C = input[3];
  params.Norg_FOM = input[4];

  return(params);
}


// [[Rcpp::export]]
double uptake_N(double N,   // UNITS: C kg
                double T,       // UNITS: 'C
                double N_limit,
                double k,
                double SWC,     // UNITS: %
                double SWC_sat) {

  // Concentration
  double u = std::max(k * pow(N, 8) / (pow(N_limit, 8) + pow(N, 8)), 0.0);
  // Temperature
  double u_t = T+20/55;
  // Water
  double u_w = std::max(pow(SWC, 8) / (pow(SWC_sat, 8) + pow(SWC, 8)), 0.0);
  double all = u * u_t * u_w;
  return(all);
}


// [[Rcpp::export]]
double uptake_C(double C,     // UNITS: C kg
                double T,     // UNITS: 'C
                double C_limit,
                double k,
                double SWC,   // UNITS: %
                double SWC_k) {

  double u = k * pow(C, 8) / (pow(C_limit, 8) + pow(C, 8)) - k;
  // Temperature
  double u_t = (T+20)/55;
  // Water
  double u_w = SWC_k * pow(SWC, 8) / (pow(0.5, 8) + pow(SWC, 8));
  double all = u * u_t * u_w;
  return(all);
}

// [[Rcpp::export]]
Rcpp::List Plant_N_Uptake(double T,
                          double SWC,
                          double m,
                          double NH4_in,
                          double NO3_in,
                          double FOM_in,
                          std::vector<double> N_limits_R,
                          std::vector<double> N_k_R,
                          std::vector<double> SWC_k_R,
                          std::vector<double> parameters,
                          double demand) {

  // Input the parameters!
  N_balence N_limits = vector_to_N_balence(N_limits_R);
  N_balence N_k = vector_to_N_balence(N_k_R);
  N_balence SWC_k = vector_to_N_balence(SWC_k_R);

  // STEP 1: Pure uptake
  // Assume that the organic and inorganic uptake is parallel
  // Inorganic NH4 affects NO3, in this function as it is a specific plant effect
  double NH4_effect_on_NO3 = parameters[0] * pow(NH4_in, 8) / (pow(parameters[1], 8) + pow(NH4_in, 8));

  // All possible N to root with NH4 modifier for NO3
  double N_to_root = (uptake_N(FOM_in, T, N_limits.Norg, N_k.Norg, SWC, SWC_k.Norg) +
                      uptake_N(NH4_in, T, N_limits.NH4, N_k.NH4, SWC, SWC_k.NH4) +
                      NH4_effect_on_NO3*uptake_N(NO3_in, T, N_limits.NO3, N_k.NO3, SWC, SWC_k.NO3)) * demand;

  return(Rcpp::List::create(Rcpp::_["N_to_plant"] = N_to_root,
                            Rcpp::_["NH4_used"] = uptake_N(NH4_in, T, N_limits.NH4, N_k.NH4, SWC, SWC_k.NH4) * demand,
                            Rcpp::_["NO3_used"] = NH4_effect_on_NO3 * uptake_N(NO3_in, T, N_limits.NO3, N_k.NO3, SWC, SWC_k.NO3) * demand,
                            Rcpp::_["Norg_used"] = uptake_N(FOM_in, T, N_limits.Norg, N_k.Norg, SWC, SWC_k.Norg) * demand));
}

// [[Rcpp::export]]
Rcpp::List Fungal_N_Uptake(double T,
                           double SWC,
                           double NH4,
                           double NO3,
                           double FOM_Norg,
                           std::vector<double> N_limits_R,
                           std::vector<double> N_k_R,
                           std::vector<double> SWC_k_R,
                           double demand) {

  N_balence N_limits = vector_to_N_balence(N_limits_R);
  N_balence N_k = vector_to_N_balence(N_k_R);
  N_balence SWC_k = vector_to_N_balence(SWC_k_R);

  double N_fungal_uptake = (uptake_N(FOM_Norg, T, N_limits.Norg, N_k.Norg, SWC, SWC_k.Norg) +
                            uptake_N(NH4, T, N_limits.NH4, N_k.NH4, SWC, SWC_k.NH4) +
                            uptake_N(NO3, T, N_limits.NO3, N_k.NO3, SWC, SWC_k.NO3))*demand;

  // TODO: consider renaming this uptake rather than used
  return Rcpp::List::create(Rcpp::_["N_to_fungal"] = N_fungal_uptake,
                            Rcpp::_["NH4_used"] = uptake_N(NH4, T, N_limits.NH4, N_k.NH4, SWC, SWC_k.NH4)*demand,
                            Rcpp::_["NO3_used"] = uptake_N(NO3, T, N_limits.NO3, N_k.NO3, SWC, SWC_k.NO3)*demand,
                            Rcpp::_["Norg_used"] = uptake_N(FOM_Norg, T, N_limits.Norg, N_k.Norg, SWC, SWC_k.Norg)*demand);
}


// [[Rcpp::export]]
Rcpp::List Microbe_Uptake(double C_microbe,                   // UNITS: C kg
                          double N_micorbe,                   // UNITS: C kg eq
                          double C_soil_compartment,
                          double NC_microbe_opt,              // UNITS: %
                          double NH4_avaliable,               // UNITS: C kg eq
                          double NO3_avaliable,               // UNITS: C kg eq
                          double Norg_avaliable,              // UNITS: C kg eq
                          double T,                           // UNITS: 'C
                          double SWC,                         // UNITS: %
                          double NC_Litter,
                          double imobilisation,
                          double assimilation,
                          std::vector<double> N_limits_R,
                          std::vector<double> N_k_R,
                          std::vector<double> SWC_k_R,
                          bool SOM_decomposers,
                          std::vector<double> respiration_microbes_params) {


  /*
   * Nitrogen limitation
   *
   * Mass limitation is in the soil model!
   */

  N_balence N_limits = vector_to_N_balence(N_limits_R);
  N_balence N_k = vector_to_N_balence(N_k_R);
  N_balence SWC_k = vector_to_N_balence(SWC_k_R);

  double Norg_uptake = uptake_N(Norg_avaliable, T, N_limits.Norg, N_k.Norg, SWC, SWC_k.Norg);   // UNITS: C kg
  double NH4_uptake = uptake_N(NH4_avaliable, T, N_limits.NH4, N_k.NH4, SWC, SWC_k.NH4);        // UNITS: C kg
  double NO3_uptake = uptake_N(NO3_avaliable, T, N_limits.NO3, N_k.NO3, SWC, SWC_k.NO3);        // UNITS: C kg

  double carbon_limitation = Norg_uptake;
  // TODO: add respiration when I have the parameters after Christmas
  // 0.2 = respiration(T, respiration_microbes_params[1], respiration_microbes_params[2]), which is the old respiration function
  // Update the parameters in the CASSIA respiration function!
  double nitrogen_limitation = (imobilisation*(NH4_uptake + NO3_uptake) + NC_microbe_opt*0.2*C_microbe)/
    (NC_Litter - NC_microbe_opt);

  /*

  double NH4_uptaken;     // UNITS: C kg
  double NO3_uptaken;     // UNITS: C kg
  double Norg_uptaken;    // UNITS: C kg
  double C_uptaken;       // UNITS: C kg

  double total_N_uptake = std::min(carbon_limitation, nitrogen_limitation);

  if (total_N_uptake == carbon_limitation) {
    C_uptaken = Norg_uptake * NC_Litter;
    Norg_uptaken = Norg_uptake;
    NH4_uptaken = NH4_uptake * (carbon_limitation / nitrogen_limitation);
    NO3_uptaken = NO3_uptake * (carbon_limitation / nitrogen_limitation);
  } else {
    C_uptaken = Norg_uptake * NC_Litter * (nitrogen_limitation / carbon_limitation);
    Norg_uptaken = Norg_uptake * (nitrogen_limitation / carbon_limitation);
    NH4_uptaken = NH4_uptake;
    NO3_uptaken = NO3_uptake;
  }

  total_N_uptake = std::max(total_N_uptake, 0.0);

  double total_N_uptaken = (NC_microbe_opt*respiration(T, respiration_microbes_params[1], respiration_microbes_params[2])*C_microbe + NC_Litter - NC_microbe_opt) * total_N_uptake;
  double Extra_FOM_uptaken = 0;
  if (SOM_decomposers) {
    total_N_uptaken =  total_N_uptaken + assimilation * C_microbe;
    // TODO fom condition?
    Extra_FOM_uptaken = 1;
  }
  */

  return(Rcpp::List::create(Rcpp::_["NH4_uptaken"] = 0, // NH4_uptaken,                                  // UNITS: C kg
                            Rcpp::_["NO3_uptaken"] = 0, // NO3_uptaken,                                  // UNITS: C kg
                            Rcpp::_["Norg_uptaken"] = 0, // Norg_uptaken,                                // UNITS: C kg
                            Rcpp::_["C_uptaken"] = 0, // C_uptaken,
                            Rcpp::_["Extra_FOM_uptaken"] = 0)); // Extra_FOM_uptaken));                                    // UNITS: C kg
}





