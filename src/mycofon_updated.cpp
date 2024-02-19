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


MYCOFON_function_out mycofon_balence(double C_biomass,
                                     double C_roots,
                                     double C_fungal,
                                     double C_ecto,
                                     double C_roots_NonStruct,
                                     double N_roots_NonStruct,
                                     double C_fungal_NonStruct,
                                     double N_fungal_NonStruct,
                                     double max_C_from_CASSIA,
                                     parameters_soil parameters_in,
                                     double NH4,
                                     double NO3,
                                     double FOM_Norg,
                                     double T,
                                     double Tmb,
                                     double SWC,
                                     bool mycofon_stratergy) {


  // Mycorrhizal rate
  // TODO: this is now in the wrong place I think!
  double m = 0.9; // C_fungal / (C_roots * parameters_in.optimal_root_fungal_biomass_ratio);
  if (m > 1) {
    std::cout << "Warning: The mycorhization value (m) is more than 1. m = " << m << " Check the fungal_mass, root_mass and optimal_root_fungal_biomass_ratio values\n";
    m = 1;
  } else if (m < 0) {
    std::cout << "Warning: The mycorhization value (m) is less than 0. m = " << m << " Check the fungal_mass, root_mass and optimal_root_fungal_biomass_ratio values\n";
    m = 0;
  }

  /*
   * Biomass accumulation
   */

  C_biomass = C_biomass + C_roots - m * C_biomass * 1/(365*2) - (1-m) * C_biomass * 1/365; // TODO: need to un hard code the turnover
  // Respiration is in CASSIA

  /*
   * DECISION VALUES CHOSEN FROM THE LAST INTERATION DATA!
   *
   * Fungal logic considered first so the N allocation value can be directly used in the plant decision function - this means that the code looks nicer!
   */

  Rcpp::List myco_decision_out = myco_decision(N_fungal_NonStruct,
                                               C_roots_NonStruct,
                                               N_roots_NonStruct,
                                               parameters_in.NC_fungal_opt);

  double N_demand_mycofon = myco_decision_out[0];
  double N_given_mycofon = myco_decision_out[1];
  double N_demand_franklin = myco_decision_out[2];
  double N_given_franklin = myco_decision_out[3];

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
  Rcpp::List plant_decision_out = plant_decision(C_roots_NonStruct,
                                                 N_roots_NonStruct,
                                                 C_fungal_NonStruct,
                                                 parameters_in.optimal_root_fungal_biomass_ratio,
                                                 m);
  // TODO: Think about this more! The last argument is the maximum value that is transferred from CASSIA
  // TODO: should this have biomass?
  double C_demand_mycofon = plant_decision_out[0];
  double C_given_mycofon = plant_decision_out[1];
  double C_demand_franklin = plant_decision_out[2];
  double C_given_franklin = plant_decision_out[3];
  double C_exudes_mycofon = plant_decision_out[4];
  double C_exudes_franklin = plant_decision_out[5];


  double C_given, plant_demand, C_exudes_plant;
  if (mycofon_stratergy) {
    C_given = C_given_mycofon;
    plant_demand = C_demand_mycofon;
    C_exudes_plant = C_exudes_mycofon;
  } else {
    C_given = C_given_franklin;
    plant_demand = C_demand_franklin;
    C_exudes_plant = C_exudes_franklin;
  }

  /*
   * UPTAKE VALUES CALCULATED
   *
   * This should be after the demand value, in case the demand term is changed, but at the moment it is always 1 so there is not a big difference!
   */

  //dN^f/dt

  Rcpp::List uptake_fungal = Fungal_N_Uptake(Tmb,
                                             SWC,
                                             NH4,
                                             NO3,
                                             FOM_Norg,
                                             parameters_in.N_limits_myco,
                                             parameters_in.N_k_myco,
                                             parameters_in.SWC_limits_myco,
                                             fungal_demand);
  double uptake_fungal_all = uptake_fungal[0];
  uptake_fungal_all = uptake_fungal_all * C_fungal;
  double uptake_fungal_NH4 = uptake_fungal[1];
  uptake_fungal_NH4 = uptake_fungal_NH4 * C_fungal;
  double uptake_fungal_NO3 = uptake_fungal[2];
  uptake_fungal_NO3 = uptake_fungal_NO3 * C_fungal;
  double uptake_fungal_Norg = uptake_fungal[3];
  uptake_fungal_Norg = uptake_fungal_Norg * C_fungal;

  // TODO: now the NH4, NO3 and FOM_Norg has changed, surely the pools need to be updated

  //dN^r/dt
  Rcpp::List uptake_plant = Plant_N_Uptake(Tmb,
                                           SWC,
                                           m,
                                           NH4,
                                           NO3,
                                           FOM_Norg,
                                           parameters_in.N_limits_plant,
                                           parameters_in.N_k_plant,
                                           parameters_in.SWC_limits_plant,
                                           parameters_in.NH4_on_NO3,
                                           plant_demand);

  double uptake_plant_all = uptake_plant[0];
  uptake_plant_all = uptake_plant_all * C_biomass;
  double uptake_plant_NH4 = uptake_plant[1];
  uptake_plant_NH4 = uptake_plant_NH4 * C_biomass;
  double uptake_plant_NO3 = uptake_plant[2];
  uptake_plant_NO3 = uptake_plant_NO3 * C_biomass;
  double uptake_plant_Norg = uptake_plant[3];
  uptake_plant_Norg = uptake_plant_Norg * C_biomass;

  /*
   * BALANCES OF THE STATE VARIABLES
   */

  double to_CASSIA = 0.02*N_roots_NonStruct; // TODO: this should link to the litter model!

  /*
   * Biomass mycorrhiza
   */

  // CALCULATE THE PARAMETERS
  // Nitrogen!

  double myco_growth_C = myco_growth(C_fungal_NonStruct, N_fungal_NonStruct, C_fungal, C_ecto, parameters_in.growth_C, parameters_in.growth_N, parameters_in.NC_fungal_opt)[1];
  double myco_growth_N = myco_growth(C_fungal_NonStruct, N_fungal_NonStruct, C_fungal, C_ecto, parameters_in.growth_C, parameters_in.growth_N, parameters_in.NC_fungal_opt)[2];

  C_fungal = C_fungal +
    myco_growth_C -
    parameters_in.turnover_mantle*C_fungal*0.5 - parameters_in.turnover_ERM*C_fungal*0.5 -
    0.00958*C_fungal; // TODO: should have a better logic for respiration

  /*
   * Non structural elements!
   */
  // TODO: need an exude percentage!
  double exudes_fungal = 0.2*C_fungal_NonStruct;
  C_fungal_NonStruct = C_fungal_NonStruct + C_given - myco_growth_C - exudes_fungal;

  N_fungal_NonStruct = N_fungal_NonStruct + uptake_fungal_all - myco_growth_N - N_given;

  N_roots_NonStruct = N_roots_NonStruct + N_given + uptake_plant_all - to_CASSIA;

  C_roots_NonStruct = C_roots_NonStruct - C_exudes_plant;

   // the to_CASSIA term should be thought to include both the nitrogen that goes into growth
   // as well as the non structural nitrogen that goes to the rest of the plant

  /*
   * OUTPUT!
   */

  MYCOFON_function_out out;
  out.C_biomass = C_biomass;
  out.C_roots = C_roots;
  out.C_fungal = C_fungal;
  out.C_roots_NonStruct = C_roots_NonStruct; // TODO: change this output I think it results in the constant values!
  out.C_fungal_NonStruct = C_fungal_NonStruct;
  out.N_roots_NonStruct = N_roots_NonStruct;
  out.N_fungal_NonStruct = N_fungal_NonStruct;
  out.uptake_plant = uptake_plant_all;
  out.uptake_NH4_plant = uptake_plant_NH4;
  out.uptake_NO3_plant = uptake_plant_NO3;
  out.uptake_Norg_plant = uptake_plant_Norg;
  out.uptake_fungal = uptake_fungal_all;
  out.uptake_NH4_fungal = uptake_fungal_NH4;
  out.uptake_NO3_fungal = uptake_fungal_NO3;
  out.uptake_Norg_fungal = uptake_fungal_Norg;
  out.from_CASSIA = max_C_from_CASSIA;
  out.to_CASSIA = to_CASSIA;
  out.Plant_demand = plant_demand;
  out.Fungal_demand = fungal_demand;
  out.Plant_given = C_given;
  out.Fungal_given = N_given;
  out.exudes_plant = C_exudes_plant;
  out.exudes_fungal = exudes_fungal;
  return(out);

}

/*
 * Mycofon Balence
 */

// [[Rcpp::export]]
Rcpp::List mycofon_balence(double C_biomass,
                           double C_roots,
                           double C_fungal,
                           double C_ecto,
                           double C_roots_NonStruct,
                           double N_roots_NonStruct,
                           double C_fungal_NonStruct,
                           double N_fungal_NonStruct,
                           double max_C_from_CASSIA,
                           std::vector<double> parameters_R,
                           double NH4,
                           double NO3,
                           double FOM_Norg,
                           double T,
                           double Tmb,
                           double SWC,
                           bool mycofon_stratergy) {

  parameters_soil parameters_in = parameters_initalise_test(parameters_R);

  MYCOFON_function_out MYCOFON_output = mycofon_balence(C_biomass,
                                                        C_roots,
                                                        C_fungal,
                                                        C_ecto,
                                                        C_roots_NonStruct,
                                                        N_roots_NonStruct,
                                                        C_fungal_NonStruct,
                                                        N_fungal_NonStruct,
                                                        max_C_from_CASSIA,
                                                        parameters_in,
                                                        NH4,
                                                        NO3,
                                                        FOM_Norg,
                                                        T,
                                                        Tmb,
                                                        SWC,
                                                        mycofon_stratergy);

  Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::_["C_biomass"] = MYCOFON_output.C_biomass,
                                               Rcpp::_["C_roots"] = MYCOFON_output.C_roots,
                                               Rcpp::_["C_fungal"] = MYCOFON_output.C_fungal,
                                               Rcpp::_["C_roots_NonStruct"] = MYCOFON_output.C_roots_NonStruct,
                                               Rcpp::_["C_fungal_NonStruct"] = MYCOFON_output.C_fungal_NonStruct,
                                               Rcpp::_["N_roots_NonStruct"] = MYCOFON_output.N_roots_NonStruct,
                                               Rcpp::_["N_fungal_NonStruct"] = MYCOFON_output.N_fungal_NonStruct,
                                               Rcpp::_["uptake_plant"] = MYCOFON_output.uptake_plant,
                                               Rcpp::_["uptake_NH4_plant"] = MYCOFON_output.uptake_NH4_plant,
                                               Rcpp::_["uptake_NO3_plant"] = MYCOFON_output.uptake_NO3_plant,
                                               Rcpp::_["uptake_Norg_plant"] = MYCOFON_output.uptake_Norg_plant,
                                               Rcpp::_["uptake_fungal"] = MYCOFON_output.uptake_fungal,
                                               Rcpp::_["uptake_NH4_fungal"] = MYCOFON_output.uptake_NH4_fungal,
                                               Rcpp::_["uptake_NO3_fungal"] = MYCOFON_output.uptake_NO3_fungal,
                                               Rcpp::_["uptake_Norg_fungal"] = MYCOFON_output.uptake_Norg_fungal,
                                               Rcpp::_["from_CASSIA"] = MYCOFON_output.from_CASSIA,
                                               Rcpp::_["to_CASSIA"] = MYCOFON_output.to_CASSIA);
  return(df);
}






