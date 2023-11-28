#include "mycomodel.h"

respiration_parameters respiration_vector_to_struct(std::vector<double> input) {
  respiration_parameters out;
  out.plant_a = input[0];
  out.plant_b = input[1];
  out.fungal_a = input[2];
  out.fungal_b = input[3];
  out.micorbe_a = input[4];
  out.micorbe_b = input[5];
  return(out);
}

// [[Rcpp::export]]
Rcpp::List mycofon_balence(double C_roots,
                           double N_roots,
                           double optimal_root_fungal_biomass_ratio,
                           double C_fungal,
                           double N_fungal,
                           double turnover_roots,
                           double turnover_roots_mycorrhized,
                           double turnover_mantle,
                           double turnover_ERM,
                           std::vector<double> respiration_parameters_R,
                           double NH4,
                           double NO3,
                           double FOM_Norg,
                           double NC_in_fungal_opt,
                           double T,
                           double Tsb,
                           double SWC,
                           std::vector<double> N_limits_Plant,
                           std::vector<double> N_k_Plant,
                           std::vector<double> SWC_k_Plant,
                           std::vector<double> N_limits_Fungal,
                           std::vector<double> N_k_Fungal,
                           std::vector<double> SWC_k_Fungal,
                           double mantle_mass,
                           double ERM_mass,
                           std::vector<double> parameters_NH4_on_NO3,
                           double growth_C,
                           double growth_N,
                           double max_C_allocation_CASSIA,
                           double allocation_N_to_rest_of_plant,
                           bool mycofon_stratergy) {
  
  /*
   * Initialise parameters
   */
  respiration_parameters respiration_params = respiration_vector_to_struct(respiration_parameters_R);
  
  // Mycorrhizal rate
  double m = C_fungal / (C_roots * optimal_root_fungal_biomass_ratio);
  if (m > 1) {
    m = 1;
    std::cout << "Warning: The mycorhization value (m) is more than 1. Check the fungal_mass, root_mass and optimal_root_fungal_biomass_ratio values\n"; 
  } else if (m < 0) {
    m = 0;
    std::cout << "Warning: The mycorhization value (m) is less than 0. Check the fungal_mass, root_mass and optimal_root_fungal_biomass_ratio values\n"; 
  }
  
  /*
   * DECISION VALUES CHOSEN FROM THE LAST INTERATION DATA!
   * 
   * Fungal logic considered first so the N allocation value can be directly used in the plant decision function - this means that the code looks nicer!
   */
  
  double N_given_mycofon = myco_decision(C_fungal,
                                         N_fungal,
                                         C_roots,
                                         N_roots,
                                         NC_in_fungal_opt,
                                         growth_C,
                                         growth_N)[1];
  
  double N_demand_mycofon = myco_decision(C_fungal,
                                         N_fungal,
                                         C_roots,
                                         N_roots,
                                         NC_in_fungal_opt,
                                         growth_C,
                                         growth_N)[0];
  
  double N_given_franklin = myco_decision(C_fungal,
                                          N_fungal,
                                          C_roots,
                                          N_roots,
                                          NC_in_fungal_opt,
                                          growth_C,
                                          growth_N)[3];
  
  double N_demand_franklin = myco_decision(C_fungal,
                                          N_fungal,
                                          C_roots,
                                          N_roots,
                                          NC_in_fungal_opt,
                                          growth_C,
                                          growth_N)[2];
  double N_given;
  double fungal_demand;
  if (mycofon_stratergy) {
    N_given = N_given_mycofon;
    fungal_demand = N_demand_mycofon;
  } else {
    N_given = N_given_franklin;
    fungal_demand = N_demand_franklin;
  }
  
  double C_given_mycofon = plant_decision(C_roots,
                                          N_roots,
                                          C_fungal,
                                          optimal_root_fungal_biomass_ratio,
                                          N_given,
                                          max_C_allocation_CASSIA)[1];
  
  double C_demand_mycofon = plant_decision(C_roots,
                                          N_roots,
                                          C_fungal,
                                          optimal_root_fungal_biomass_ratio,
                                          N_given,
                                          max_C_allocation_CASSIA)[0];
  
  double C_given_franklin = plant_decision(C_roots,
                                           N_roots,
                                           C_fungal,
                                           optimal_root_fungal_biomass_ratio,
                                           N_given,
                                           max_C_allocation_CASSIA)[3];
  
  double C_demand_franklin = plant_decision(C_roots,
                                           N_roots,
                                           C_fungal,
                                           optimal_root_fungal_biomass_ratio,
                                           N_given,
                                           max_C_allocation_CASSIA)[2];
  
  double C_given;
  double plant_demand; 
  if (mycofon_stratergy) {
    C_given = C_given_mycofon;
    plant_demand = C_demand_mycofon;
  } else {
    C_given = C_given_franklin;
    plant_demand = C_demand_franklin;
  }
  
  
  /*
   * UPTAKE VALUES CALCULATED
   * 
   * This should be after the demand value, in case the demand term is changed, but at the moment it is always 1 so there is not a big difference!
   */
  
  //dN^f/dt
  double uptake_fungal = Fungal_N_Uptake(T,
                                         SWC,
                                         NH4,
                                         NO3,
                                         FOM_Norg,
                                         N_limits_Fungal,
                                         N_k_Fungal,
                                         SWC_k_Fungal,
                                         fungal_demand)[0];
  
  //dN^r/dt
  double uptake_plant = Plant_N_Uptake(T,
                                       SWC,
                                       m, 
                                       NH4,
                                       NO3,
                                       FOM_Norg,
                                       N_limits_Plant,
                                       N_k_Plant,
                                       SWC_k_Plant,
                                       parameters_NH4_on_NO3,
                                       plant_demand)[0];
  
  /*
   * BALANCES OF THE STATE VARIABLES
   */
  
  // CALCULATE THE PARAMETERS
  
  // Nitrogen!
  double root_NC_ratio = N_roots / C_roots;
  double fungal_NC_ratio = N_fungal / C_fungal;
  
  double from_CASSIA = max_C_allocation_CASSIA; // could change this to be a function that would work with the gradient that is happening in CASSIA
  double to_CASSIA = allocation_N_to_rest_of_plant; // TODO: sort this!
  
  double myco_growth_C = myco_growth(C_fungal, N_fungal, growth_C, growth_N)[1]; // TODO: currently C_roots rather than sugar as the model isn't connected in that way
  double myco_growth_N = myco_growth(C_fungal, N_fungal, growth_C, growth_N)[2];
  
  N_roots = N_roots +
    N_given + 
    uptake_plant*C_roots - // This should be by the biomass, so decided that it is multiplied by C^r rather than N^r, maybe should actually be a surface area equation
    (1 - m)*C_roots*turnover_roots*root_NC_ratio - 
    m*C_roots*turnover_roots_mycorrhized*root_NC_ratio -
    to_CASSIA;
  
  N_fungal = N_fungal + 
    uptake_fungal*C_fungal - 
    myco_growth_N -
    (0.5*turnover_mantle + 0.5*turnover_ERM)*C_fungal*fungal_NC_ratio - 
    N_given;
  
  // dC^r/dt TODO: this need to be linked with the CASSIA C sections, also link the amount going to CASSIA
  C_roots = C_roots + 
    from_CASSIA - 
    (1 - m)*C_roots*turnover_roots - 
    m*C_roots*turnover_roots_mycorrhized - 
    respiration(Tsb, respiration_params.plant_a, respiration_params.plant_b)*C_roots - 
    C_given;

  // dC^f/dt
  // TODO: what is the percentage_C_biomass doing?
  C_fungal = C_fungal + 
    C_given - 
    myco_growth_C -
    turnover_mantle*mantle_mass - turnover_ERM*ERM_mass -
    respiration(Tsb, respiration_params.fungal_a, respiration_params.fungal_b)*C_fungal;
  
  /*
   * OUTPUT!
   */
  
  return Rcpp::List::create(Rcpp::_["C_roots"] = C_roots,
                            Rcpp::_["C_fungal"] = C_fungal,
                            Rcpp::_["N_roots"] = N_roots,
                            Rcpp::_["N_fungal"] = N_fungal,
                            Rcpp::_["uptake_plant"] = uptake_plant*C_roots,
                            Rcpp::_["uptake_fungal"] = uptake_fungal*N_roots,
                            Rcpp::_["from_CASSIA"] = from_CASSIA,
                            Rcpp::_["to_CASSIA"] = to_CASSIA);

}


