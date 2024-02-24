#include "CASSIA.h"


// Input R values in the correct structure - N_balence
N_balence vector_to_N_balence(std::vector<double> input) {
  N_balence params;

  params.NH4 = input[0];
  params.NO3 = input[1];
  params.Norg = input[2];

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
                double SWC,     // UNITS: %
                double N_limit,
                double k,
                double SWC_limit) {

  double u, u_t, u_w, all;
  if (N < 0) {
    std::cerr << "Warning! Pool is less than 0.";
    all = 0;
  } else {
    // Concentration
    u = std::max(k * pow(N/N_limit, 8.0) / (1 + pow(N/N_limit, 8.0)), 0.0);
    // Temperature
    u_t = std::max(T/30.0, 0.0); // TODO: these parameters shouldn't be hard coded
    // Water
    u_w = std::max(pow(SWC, 8.0) / (pow(SWC_limit, 8.0) + pow(SWC, 8.0)), 0.0);
    all = u * u_t * u_w;
    if (u_t < 0 | u_t > 1) {
      std::cout << "Warning! u_t out of bounds. u_t = " << u_t << "\n";
    }
    if (u_w < 0 | u_w > 1) {
      std::cout << "Warning! u_t out of bounds. u_w = " << u_w << "\n";
    }
  }

  if (N < all) {
    std::cerr << "Warning! Uptake is larger than the N pool: set to the N pool ";
    all = N;
  }
  if (all != all) {
    // std::cerr << " Uptake gives NA: set to 0. ";
    all = 0.0;
  }
  if (all < 0.0) {
    std::cerr << "Warning! Uptake is less than 0! \n";
    all = 0.0;
  }
  return(all);
}


// [[Rcpp::export]]
double uptake_C(double C,     // UNITS: C kg
                double T,     // UNITS: 'C
                double SWC,   // UNITS: %
                double C_limit,
                double k,
                double SWC_limit) {

  double u, u_t, u_w, all;
  if (C < 0) {
    std::cerr << "Warning! Pool is less than 0.";
    all = 0;
  } else {
    // Concentration
    u = std::max(k * pow(C/C_limit, 8.0) / (1 + pow(C/C_limit, 8.0)), 0.0);
    // Temperature
    u_t = std::max(T/30.0, 0.0);
    // Water
    u_w = pow(SWC, 8.0) / (pow(SWC_limit, 8.0) + pow(SWC, 8.0));
    all = u * u_t * u_w;
    if (u_t < 0 | u_t > 1) {
      std::cout << "Warning! u_t out of bounds. u_t = " << u_t << "\n";
    }
    if (u_w < 0 | u_w > 1) {
      std::cout << "Warning! u_t out of bounds. u_w = " << u_w << "\n";
    }
  }

  if (C < all) {
    std::cerr << "Warning! Uptake is larger than the N pool: set to the N pool ";
    all = C;
  }
  if (all != all) {
    // std::cerr << " Uptake gives NA: set to 0. Check the weather input!\n";
    all = 0.0;
  }
  if (all < 0.0) {
    std::cerr << "Warning! Uptake is less than 0! \n";
    all = 0.0;
  }
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
                          std::vector<double> SWC_limits_R,
                          double NH4_on_NO3,
                          double demand) {

  // Input the parameters!
  N_balence N_limits = vector_to_N_balence(N_limits_R);
  N_balence N_k = vector_to_N_balence(N_k_R);
  N_balence SWC_limits = vector_to_N_balence(SWC_limits_R);

  // STEP 1: Pure uptake
  // Assume that the organic and inorganic uptake is parallel
  // Inorganic NH4 affects NO3, in this function as it is a specific plant effect
  double NH4_effect_on_NO3 = pow(NH4_in/NH4_on_NO3, 8.0) / (1.0 + pow(NH4_in/NH4_on_NO3, 8.0));

  // All possible N to root with NH4 modifier for NO3
  double N_to_root = (1 - m) * (uptake_N(FOM_in, T, SWC,N_limits.Norg, N_k.Norg, SWC_limits.Norg) +
                      uptake_N(NH4_in, T, SWC, N_limits.NH4, N_k.NH4, SWC_limits.NH4) +
                      NH4_effect_on_NO3*uptake_N(NO3_in, T, SWC, N_limits.NO3, N_k.NO3, SWC_limits.NO3)) * demand;

  return(Rcpp::List::create(Rcpp::_["N_to_plant"] = N_to_root,
                            Rcpp::_["NH4_used"] = (1 - m) * uptake_N(NH4_in, T, SWC, N_limits.NH4, N_k.NH4, SWC_limits.NH4) * demand,
                            Rcpp::_["NO3_used"] = (1 - m) * NH4_effect_on_NO3 * uptake_N(NO3_in, T, SWC, N_limits.NO3, N_k.NO3, SWC_limits.NO3) * demand,
                            Rcpp::_["Norg_used"] = (1 - m) * uptake_N(FOM_in, T, SWC, N_limits.Norg, N_k.Norg, SWC_limits.Norg) * demand));
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

  double N_fungal_uptake = (uptake_N(FOM_Norg, T, SWC, N_limits.Norg, N_k.Norg, SWC_k.Norg) +
                            uptake_N(NH4, T, SWC, N_limits.NH4, N_k.NH4, SWC_k.NH4) +
                            uptake_N(NO3, T, SWC, N_limits.NO3, N_k.NO3, SWC_k.NO3))*demand;

  // TODO: consider renaming this uptake rather than used
  return Rcpp::List::create(Rcpp::_["N_to_fungal"] = N_fungal_uptake,
                            Rcpp::_["NH4_used"] = uptake_N(NH4, T, SWC, N_limits.NH4, N_k.NH4, SWC_k.NH4)*demand,
                            Rcpp::_["NO3_used"] = uptake_N(NO3, T, SWC, N_limits.NO3, N_k.NO3, SWC_k.NO3)*demand,
                            Rcpp::_["Norg_used"] = uptake_N(FOM_Norg, T, SWC, N_limits.Norg, N_k.Norg, SWC_k.Norg)*demand);
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
                          double FOM_Norg) {


  /*
   * Nitrogen limitation
   *
   * Mass limitation is in the soil model!
   */

  N_balence N_limits = vector_to_N_balence(N_limits_R);
  N_balence N_k = vector_to_N_balence(N_k_R);
  N_balence SWC_k = vector_to_N_balence(SWC_k_R);

  double NH4_uptake = uptake_N(NH4_avaliable, T, SWC, N_limits.NH4, N_k.NH4, SWC_k.NH4);        // UNITS: C kg
  double NO3_uptake = uptake_N(NO3_avaliable, T, SWC, N_limits.NO3, N_k.NO3, SWC_k.NO3);        // UNITS: C kg
  double Norg_uptake = uptake_N(Norg_avaliable, T, SWC, N_limits.Norg, N_k.Norg, SWC_k.Norg);   // UNITS: C kg

  double carbon_limitation = C_microbe * Norg_uptake;
  // TODO: write the equations in the other document
  // TODO: the decomposition is assumed to release exactly the same amount of nitrogen and carbon currently!

  // TODO: add respiration when I have the parameters after Christmas
  // 0.2 = respiration(T, respiration_microbes_params[1], respiration_microbes_params[2]), which is the old respiration function
  // Update the parameters in the CASSIA respiration function!
  double nitrogen_limitation = (imobilisation * (NH4_uptake + NO3_uptake) + NC_microbe_opt * 0.2 * C_microbe)/
    (NC_Litter - NC_microbe_opt);

  double NH4_uptaken;     // UNITS: C kg
  double NO3_uptaken;     // UNITS: C kg
  double Norg_uptaken;    // UNITS: C kg
  double C_uptaken;       // UNITS: C kg

  double total_N_uptake = std::max(0.0, std::min(carbon_limitation, nitrogen_limitation));
  // NOTE: this only works because the carbon and nitrogen have the same units

  if (total_N_uptake == carbon_limitation) {
    C_uptaken = Norg_uptake * NC_Litter;
    Norg_uptaken = Norg_uptake;
    NH4_uptaken = NH4_uptake * (carbon_limitation / nitrogen_limitation);
    NO3_uptaken = NO3_uptake * (carbon_limitation / nitrogen_limitation);
  } else if (total_N_uptake == nitrogen_limitation) {
    C_uptaken = Norg_uptake * NC_Litter * (nitrogen_limitation / carbon_limitation);
    Norg_uptaken = Norg_uptake * (nitrogen_limitation / carbon_limitation);
    NH4_uptaken = NH4_uptake;
    NO3_uptaken = NO3_uptake;
  } else {
    C_uptaken = 0;
    Norg_uptaken = 0;
    NH4_uptaken = 0;
    NO3_uptaken = 0;
  }

  double total_N_uptaken = (NC_microbe_opt*0.2*C_microbe + NC_Litter - NC_microbe_opt) * total_N_uptake;

  double Extra_FOM_uptaken = 0;
  if (SOM_decomposers) {
    total_N_uptaken =  total_N_uptaken + assimilation * C_microbe;
    // TODO from condition?
    Extra_FOM_uptaken = uptake_N(Norg_uptake, T, SWC, N_limits.Norg, N_k.Norg, SWC_k.Norg);
  }

  return(Rcpp::List::create(Rcpp::_["NH4_uptaken"] = NH4_uptaken,               // UNITS: C kg
                            Rcpp::_["NO3_uptaken"] = NO3_uptaken,               // UNITS: C kg
                            Rcpp::_["Norg_uptaken"] = Norg_uptaken,             // UNITS: C kg
                            Rcpp::_["C_uptaken"] = C_uptaken,                   // C_uptaken
                            Rcpp::_["Extra_FOM_uptaken"] = Extra_FOM_uptaken)); // Extra_FOM_uptaken                                   // UNITS: C kg
}





