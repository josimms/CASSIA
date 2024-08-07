#include "CASSIA.h"

SYMPHONY_output symphony_multiple_FOM_daily(double Tmb,
                                            double SWC,
                                            double C_FOM_needles_old,
                                            double C_FOM_woody_old,
                                            double C_FOM_roots_old,
                                            double C_FOM_mantle_old,
                                            double C_FOM_ERM_old,
                                            double C_Exudes,
                                            double C_SOM_old,
                                            double N_SOM_old,
                                            double C_decompose_FOM,
                                            double C_decompose_SOM,
                                            double N_decompose_FOM,
                                            double N_decompose_SOM,
                                            double Litter_needles,
                                            double Litter_woody,
                                            double Litter_roots,
                                            double Litter_mantle,
                                            double Litter_ERM,
                                            double exudes_plant,
                                            double exudes_fungal,
                                            double imobilisation,
                                            double assimilation,
                                            double NH4_old,
                                            double NO3_old,
                                            double NC_needles,
                                            double NC_woody,
                                            double NC_roots,
                                            double NC_mantle,
                                            double NC_ERM,
                                            double NH4_used_Plant,
                                            double NH4_used_Fungal,
                                            double NO3_used_Plant,
                                            double NO3_used_Fungal,
                                            double FOM_Norg_used_Plant,
                                            double FOM_Norg_used_Fungal,
                                            double SOM_Norg_used,
                                            std::vector<double> N_limits_R,
                                            std::vector<double> N_k_R,
                                            std::vector<double> SWC_k_R,
                                            double NC_microbe_opt,
                                            double microbe_turnover,
                                            bool tests)
{
  /*
   * Initialization or declaration
   */


  // INITILISATION FROM VECTORS
  // Parameters into symphony_parameters format

  // SOIL INITIALISATION NEW INPUT!
  double C_FOM_needles = C_FOM_needles_old + Litter_needles;             // C kg
  double C_FOM_woody = C_FOM_woody_old + Litter_woody;                   // C kg
  double C_FOM_roots = C_FOM_roots_old + Litter_roots;                   // C kg
  double C_FOM_mantle = C_FOM_mantle_old + Litter_mantle;                // C kg
  double C_FOM_ERM = C_FOM_ERM_old + Litter_ERM;                         // C kg

  // AGGREGATION
  double C_FOM = C_FOM_needles +
    C_FOM_woody +
    C_FOM_roots +
    C_FOM_mantle +
    C_FOM_ERM; // C kg
  double N_FOM = C_FOM_needles * NC_needles +
    C_FOM_woody * NC_woody +
    C_FOM_roots * NC_roots +
    C_FOM_mantle * NC_mantle +
    C_FOM_ERM * NC_ERM; // C kg

  /*
   * Nitrogen processes
   */

  // STEP 1: Set the values, considering the other models
  double NH4 = NH4_old - NH4_used_Plant - NH4_used_Fungal;             // C kg
  double NO3 = NO3_old - NO3_used_Plant - NO3_used_Fungal;             // C kg
  N_FOM = N_FOM - FOM_Norg_used_Plant - FOM_Norg_used_Fungal;          // C kg TODO: there is a chance that this will go negative!
  double N_SOM = N_SOM_old;                                            // C kg''


  // TODO: what about the carbon considerations in the N used!

  // STEP 2: consider the uptake functions, for FOM and SOM, MICROBES
  // this includes the decomposition and mineralisation / immobilisation
  // Apart from the biomass that the microbes specialise in, the second N uptake is considered to be equally driven between the types

  // FOM
  Rcpp::List FOM_after_microbe_activity_list = Microbe_Uptake(C_decompose_FOM, N_decompose_FOM, C_Exudes, C_FOM, NC_microbe_opt, NH4, NO3, N_FOM, Tmb, SWC, imobilisation, assimilation, N_limits_R, N_k_R, SWC_k_R, false, 0, true);
  N_balence FOM_after_microbe_activity = list_to_N_balence(FOM_after_microbe_activity_list);    // C kg

  NO3 = NO3 - C_decompose_FOM*FOM_after_microbe_activity.NO3; // C kg eq
  NH4 = NH4 - C_decompose_FOM*FOM_after_microbe_activity.NH4;  // C kg eq

  N_FOM = N_FOM - C_decompose_FOM * FOM_after_microbe_activity.Norg;   // C kg eq
  N_SOM = N_SOM + microbe_turnover * (N_decompose_FOM);    // C kg eq

  // SOM
  Rcpp::List SOM_after_microbe_activity_list = Microbe_Uptake(C_decompose_SOM, N_decompose_SOM, C_Exudes, C_SOM_old, NC_microbe_opt, NH4, NO3, N_SOM, Tmb, SWC, imobilisation, assimilation, N_limits_R, N_k_R, SWC_k_R, true, N_FOM, true);
  N_balence SOM_after_microbe_activity = list_to_N_balence(SOM_after_microbe_activity_list);    // C kg

  // Update pure inorganic pools
  NO3 = NO3 - SOM_after_microbe_activity.NO3;    // C kg eq
  NH4 = NH4 - SOM_after_microbe_activity.NH4;    // C kg eq

  N_SOM = N_SOM - SOM_after_microbe_activity.Norg + microbe_turnover * (N_decompose_SOM + N_decompose_FOM);    // C kg eq
  N_FOM = N_FOM - SOM_after_microbe_activity.Norg_FOM;    // C kg eq

  /*
   * Carbon processes
   */

  // STEP 1: Update the FOM mass
  double total_decomposition = FOM_after_microbe_activity.C + SOM_after_microbe_activity.C;   // C kg
  // TODO: the nitrogen transfer should be by N and litter type rather than size of pool!
  C_FOM_needles = C_FOM_needles - (C_FOM_needles/C_FOM)*total_decomposition;            // C kg
  C_FOM_woody = C_FOM_woody - (C_FOM_woody/C_FOM)*total_decomposition;                  // C kg
  C_FOM_roots = C_FOM_roots - (C_FOM_roots/C_FOM)*total_decomposition;                  // C kg
  C_FOM_mantle = C_FOM_mantle - (C_FOM_mantle/C_FOM)*total_decomposition;               // C kg
  C_FOM_ERM = C_FOM_ERM - (C_FOM_ERM/C_FOM)*total_decomposition;                        // C kg

  C_Exudes = C_Exudes + exudes_fungal + exudes_plant -
    FOM_after_microbe_activity.C_exudes -
    SOM_after_microbe_activity.C_exudes;

  // STEP 2: Update the SOM mass
  double C_SOM = C_SOM_old + microbe_turnover * (C_decompose_FOM + C_decompose_SOM) - SOM_after_microbe_activity.C;   // C kg

  // STEP 3: Update the microbes
  double resp = 0.02; // TODO: add respiration function here!
  C_decompose_FOM = (1 - microbe_turnover - resp) * C_decompose_FOM + FOM_after_microbe_activity.C + FOM_after_microbe_activity.C_exudes;   // C kg
  N_decompose_FOM = (1 - microbe_turnover) * N_decompose_FOM + FOM_after_microbe_activity.NH4 + FOM_after_microbe_activity.NO3 + FOM_after_microbe_activity.Norg;   // C kg

  C_decompose_SOM = (1 - microbe_turnover - resp) * C_decompose_SOM + SOM_after_microbe_activity.C + SOM_after_microbe_activity.C_exudes;   // C kg
  N_decompose_SOM = (1 - microbe_turnover) * N_decompose_SOM + SOM_after_microbe_activity.NH4 + SOM_after_microbe_activity.NO3 + SOM_after_microbe_activity.Norg;   // C kg
  // TODO: 0.2 is a placeholder for respiration

  double amino_acids = 1; // TOOD: add the organic S storage into the model


  SYMPHONY_output out;
  out.C_decompose_FOM = C_decompose_FOM;          // C kg
  out.C_decompose_SOM = C_decompose_SOM;          // C kg
  out.C_FOM_ERM = C_FOM_ERM;                      // C kg
  out.C_FOM_mantle = C_FOM_mantle;                // C kg
  out.C_FOM_needles = C_FOM_needles;              // C kg
  out.C_FOM_roots = C_FOM_roots;                  // C kg
  out.C_FOM_woody = C_FOM_woody;                  // C kg
  out.C_exudes = C_Exudes;                        // C kg
  out.C_SOM = C_SOM;                              // C kg
  out.N_decompose_FOM = N_decompose_FOM;          // C kg
  out.N_decompose_SOM = N_decompose_SOM;          // C kg
  out.N_FOM = N_FOM;                              // C kg
  out.N_SOM = N_SOM;                              // C kg
  out.NC_ERM = NC_ERM;                            // C kg
  out.NC_mantle = NC_mantle;                      // C kg
  out.NC_needles = NC_needles;                    // C kg
  out.NC_roots = NC_roots;                        // C kg
  out.NC_woody = NC_woody;                        // C kg
  out.NH4 = NH4;                                  // C kg
  out.NO3 = NO3;                                  // C kg
  out.SOM_Norg_used = SOM_Norg_used;              // C kg
  out.Microbe_respiration_per_mass = resp;
  out.NH4_Uptake_Microbe_FOM = FOM_after_microbe_activity.NH4;
  out.NO3_Uptake_Microbe_FOM = FOM_after_microbe_activity.NO3;
  out.Norg_Uptake_Microbe_FOM = FOM_after_microbe_activity.Norg;
  out.C_Uptake_Microbe_FOM = FOM_after_microbe_activity.C;
  out.C_exudes_Uptake_Microbe_FOM = FOM_after_microbe_activity.C_exudes;
  out.NH4_Uptake_Microbe_SOM = SOM_after_microbe_activity.NH4;
  out.NO3_Uptake_Microbe_SOM = SOM_after_microbe_activity.NO3;
  out.Norg_Uptake_Microbe_SOM = SOM_after_microbe_activity.Norg;
  out.C_Uptake_Microbe_SOM = SOM_after_microbe_activity.C;
  out.C_exudes_Uptake_Microbe_SOM = SOM_after_microbe_activity.C_exudes;

  return(out);
}


