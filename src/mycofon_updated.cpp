#include "CASSIA.h"

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


MYCOFON_function_out mycofon_balence(double C_roots,
                                     double C_growth,
                                     double N_roots,
                                     double C_fungal,
                                     double N_fungal,
                                     parameters_soil parameters_in,
                                     double NH4,
                                     double NO3,
                                     double FOM_Norg,
                                     double T,
                                     double Tsb,
                                     double SWC,
                                     double mantle_mass,
                                     double ERM_mass,
                                     bool mycofon_stratergy) {


  // Mycorrhizal rate
  double m = 0.9; //  C_fungal / (C_roots * parameters_in.optimal_root_fungal_biomass_ratio);
  if (m > 1) {
    std::cout << "Warning: The mycorhization value (m) is more than 1. m = " << m << " Check the fungal_mass, root_mass and optimal_root_fungal_biomass_ratio values\n";
    m = 1;
  } else if (m < 0) {
    std::cout << "Warning: The mycorhization value (m) is less than 0. m = " << m << " Check the fungal_mass, root_mass and optimal_root_fungal_biomass_ratio values\n";
    m = 0;
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
                                         parameters_in.NC_fungal_opt,
                                         parameters_in.growth_C,
                                         parameters_in.growth_N)[1];
  // std::cout << "N_given_mycofon: " << N_given_mycofon;

  double N_demand_mycofon = myco_decision(C_fungal,
                                         N_fungal,
                                         C_roots,
                                         N_roots,
                                         parameters_in.NC_fungal_opt,
                                         parameters_in.growth_C,
                                         parameters_in.growth_N)[0];
  // std::cout << " N_demand_mycofon: " << N_demand_mycofon << "\n";

  double N_given_franklin = myco_decision(C_fungal,
                                          N_fungal,
                                          C_roots,
                                          N_roots,
                                          parameters_in.NC_fungal_opt,
                                          parameters_in.growth_C,
                                          parameters_in.growth_N)[3];

  double N_demand_franklin = myco_decision(C_fungal,
                                          N_fungal,
                                          C_roots,
                                          N_roots,
                                          parameters_in.NC_fungal_opt,
                                          parameters_in.growth_C,
                                          parameters_in.growth_N)[2];
  double N_given;
  double fungal_demand;
  if (mycofon_stratergy) {
    N_given = N_given_mycofon;
    fungal_demand = N_demand_mycofon;
  } else {
    N_given = N_given_franklin;
    fungal_demand = N_demand_franklin;
  }

  // NOTE! This function takes the maximum nutrient amount and then allocates the amount the plant / fungi wants to give
  double C_given_mycofon = plant_decision(C_roots,
                                          N_roots,
                                          C_fungal,
                                          parameters_in.optimal_root_fungal_biomass_ratio,
                                          N_given,
                                          0.7*C_roots)[1]; // TODO: Think about this more! The last argument is the maximum value that is transfered from CASSIA

  double C_demand_mycofon = plant_decision(C_roots,
                                          N_roots,
                                          C_fungal,
                                          parameters_in.optimal_root_fungal_biomass_ratio,
                                          N_given,
                                          0.7*C_roots)[0];

  double C_given_franklin = plant_decision(C_roots,
                                           N_roots,
                                           C_fungal,
                                           parameters_in.optimal_root_fungal_biomass_ratio,
                                           N_given,
                                           0.7*C_roots)[3];

  double C_demand_franklin = plant_decision(C_roots,
                                           N_roots,
                                           C_fungal,
                                           parameters_in.optimal_root_fungal_biomass_ratio,
                                           N_given,
                                           0.7*C_roots)[2];

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
  Rcpp::List uptake_fungal = Fungal_N_Uptake(T,
                                             SWC,
                                             NH4,
                                             NO3,
                                             FOM_Norg,
                                             parameters_in.N_limits_fungal,
                                             parameters_in.N_k_fungal,
                                             parameters_in.SWC_k_fungal,
                                             fungal_demand);
  double uptake_fungal_all = uptake_fungal[0];
  double uptake_fungal_NH4 = uptake_fungal[1];
  double uptake_fungal_NO3 = uptake_fungal[2];
  double uptake_fungal_Norg = uptake_fungal[3];

  //dN^r/dt
  Rcpp::List uptake_plant = Plant_N_Uptake(T,
                                           SWC,
                                           m,
                                           NH4,
                                           NO3,
                                           FOM_Norg,
                                           parameters_in.N_limits_plant,
                                           parameters_in.N_k_plant,
                                           parameters_in.SWC_k_plant,
                                           parameters_in.NH4_on_NO3,
                                           plant_demand);
  double uptake_plant_all = uptake_plant[0];
  double uptake_plant_NH4 = uptake_plant[1];
  double uptake_plant_NO3 = uptake_plant[2];
  double uptake_plant_Norg = uptake_plant[3];

  /*
   * BALANCES OF THE STATE VARIABLES
   */

  // CALCULATE THE PARAMETERS

  // Nitrogen!
  double root_NC_ratio = N_roots / C_roots;
  double fungal_NC_ratio = N_fungal / C_fungal;

  double myco_growth_C = myco_growth(C_fungal, N_fungal, parameters_in.growth_C, parameters_in.growth_N)[1]; // TODO: currently C_roots rather than sugar as the model isn't connected in that way
  double myco_growth_N = myco_growth(C_fungal, N_fungal, parameters_in.growth_C, parameters_in.growth_N)[2];

  double to_CASSIA = std::max(0.01*N_roots, 0.0);
  N_roots = N_roots +
    N_given +
    uptake_plant_all*C_roots - // This should be by the biomass, so decided that it is multiplied by C^r rather than N^r, maybe should actually be a surface area equation
    // TODO: Should these be in CASSIA instead?a
    //(1 - m)*C_roots*parameters_in.turnover_roots*root_NC_ratio -
    //m*C_roots*parameters_in.turnover_roots_mycorrhized*root_NC_ratio -
    to_CASSIA; // TODO: rethink: last term is the amount of N that is going from the roots to the rest of the tree

  N_fungal = N_fungal +
    uptake_fungal_all*C_fungal +
    myco_growth_N -
    //(0.5*parameters_in.turnover_mantle + 0.5*parameters_in.turnover_ERM)*C_fungal*fungal_NC_ratio - // TODO: although this is correct I need to be consistent with the fact that ERM is 50% of the fungal biomass
    N_given;

  // dC^r/dt TODO: this need to be linked with the CASSIA C sections, also link the amount going to CASSIA
  // TODO: this balance is currently biomass rather than sugar based!
  C_roots = C_roots + C_growth + 0.5*C_roots -
    C_given; // TODO: the C given should be proportional to the mycorrhized amount
    // (1 - m)*C_roots*parameters_in.turnover_roots -
    // m*C_roots*parameters_in.turnover_roots_mycorrhized -
  // 0.2 is a placeholder for respiration

  // dC^f/dt
  C_fungal = C_fungal + myco_growth_C + // TODO: should I put this in the CASSIA section?
    C_given -
    parameters_in.turnover_mantle*C_fungal*0.5 - parameters_in.turnover_ERM*C_fungal*0.5 - // TODO: should I put this in the CASSIA section?
    0.2*C_fungal;
  // 0.2 is a placeholder for respiration

  /*
   * OUTPUT!
   */

  MYCOFON_function_out out;
  out.C_roots = C_roots;
  out.C_fungal = C_fungal;
  out.N_roots = N_roots;
  out.N_fungal = N_fungal;
  out.uptake_plant = uptake_plant_all*C_roots;
  out.uptake_NH4_plant = uptake_plant_NH4*C_roots;
  out.uptake_NO3_plant = uptake_plant_NO3*C_roots;
  out.uptake_Norg_plant = uptake_plant_Norg*C_roots;
  out.uptake_fungal = uptake_fungal_all*C_fungal;
  out.uptake_NH4_fungal = uptake_fungal_NH4*C_roots;
  out.uptake_NO3_fungal = uptake_fungal_NO3*C_roots;
  out.uptake_Norg_fungal = uptake_fungal_Norg*C_roots;
  out.from_CASSIA = 0.5*C_roots;
  out.to_CASSIA = to_CASSIA;
  out.Plant_demand = plant_demand;
  out.Fungal_demand = fungal_demand;
  out.Plant_given = C_given;
  out.Fungal_given = N_given;
  return(out);

}

