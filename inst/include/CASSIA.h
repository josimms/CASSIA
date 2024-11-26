#include "CASSIA_types.h"
// CASSIA types includes function_structures.h, Rcpp.h, RcppCommon.h
// isostrem and vector are also in function_structures
#include "mycomodel.h"
// This includes soil and mycorrhizal parameters
#include <numeric>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "final_parameters.h"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <iomanip>
#include <prelesglobals.h>
#include <vector>

#ifndef CASSIA_H
#define CASSIA_H

struct Tree_out
{
  double D0;
  double h0;
  double LH0;
  double LN0;
  double LR0;
};

double emergancy(double sugar, double starch, double tau_emergancy, double lower_bound);

double storage_update(double alfa, double sugar, double starch, double Wala);

/*
 * PRELES
 */

/*
 * Yearly Initial Conditions
 */

struct yearly_in {
  // Growth
  std::vector<double> g_sH;
  std::vector<double> sH;
  std::vector<double> fH;
  std::vector<double> GH;
  std::vector<double> cumsum_GH;
  std::vector<double> HH;
  std::vector<double> height_pot_growth;
  std::vector<double> LH;
  std::vector<double> g_sN_T;
  std::vector<double> sN;
  std::vector<double> fN;
  std::vector<double> GN;
  std::vector<double> cumsum_GN;
  std::vector<double> needle_pot_growth;
  std::vector<double> LN;
  std::vector<double> HN;
  std::vector<double> g_sD_T;
  std::vector<double> sD;
  std::vector<double> fD;
  std::vector<double> GD;
  std::vector<double> diameter_pot_growth;
  std::vector<double> LD;
  std::vector<double> GPP_ref;
  std::vector<double> S_GPP;
  std::vector<double> S_GPP_ref;
  std::vector<double> dS_GPP;
  std::vector<double> dS_GPP_ref;
  std::vector<double> g_R;
  std::vector<double> sR;
  std::vector<double> fR;
  std::vector<double> GR;
  std::vector<double> roots_pot_growth;
  std::vector<double> LR;
  std::vector<double> tot_cells_pot;
  std::vector<double> tot_cells;
  std::vector<double> n_rows;
  std::vector<double> g;
  std::vector<double> sRc;

  // Xylogenesis
  std::vector<double> tau_E;
  std::vector<double> tau_W;
  std::vector<double> n_E;
  std::vector<double> n_W;
  std::vector<double> n_M;
  std::vector<double> n_E_pot;
  std::vector<double> n_W_pot;
  std::vector<double> n_M_pot;
  // std::vector<double> GD; // TOOD: this is duplicated!
  std::vector<double> CE_ew;
  std::vector<double> CE_lw;
  std::vector<double> CW;
  std::vector<double> carbon_enlargement;
  std::vector<double> wall_growth;
  std::vector<double> wall_tot;
  std::vector<double> ew_cells;
  std::vector<double> lw_cells;

  // Sugar

};

struct eternal_parameters {
  double absolute_zero; // TOOD; what should this value be?
};

yearly_in yearly_initial_conditions(double days_per_year);


/*
 * sugar_model.cpp
 */

carbo_balance sugar_model(int year,
                          int day,
                          double TAir,
                          double PF,

                          CASSIA_common common,
                          CASSIA_parameters parameters,

                          double D00,
                          double sH,
                          respiration_out resp,

                          bool sperling_sugar_model,
                          bool tree_alive,
                          bool storage_grows,
                          double needles_mass, // Repola
                          double temperature_equilibrium,

                          growth_out pot_growth,

                          carbo_tracker sugar,
                          carbo_tracker starch,

                          carbo_values_out parameters_in);

/*
 * growth.cpp
 */

growth_out growth(int day,
                  int year,
                  double TAir,
                  double TSoil_A,
                  double TSoil_B,
                  double Soil_Moisture,
                  double PF,
                  double GPP_ref,
                  bool root_as_Ding,
                  bool xylogenesis_option,
                  bool environmental_effect_xylogenesis,
                  bool sD_estim_T_count,
                  CASSIA_common common,
                  CASSIA_parameters parameters,
                  CASSIA_ratios ratio,
                  double CH,
                  double B0,
                  double GPP_mean,
                  double GPP_previous_sum,

                  bool LH_estim,
                  bool LN_estim,
                  bool LD_estim,
                  bool tests,

                  growth_values_out growth_previous,
                  double last_year_HH,
                  int no_day);

Rcpp::List growth_wrapper(int day,
                          int year,
                          double TAir,
                          double TSoil_A,
                          double TSoil_B,
                          double Soil_Moisture,
                          double PF,
                          double GPP_ref,
                          bool root_as_Ding,
                          bool xylogenesis_option,
                          bool environmental_effect_xylogenesis,
                          bool sD_estim_T_count,
                          CASSIA_common pCASSIA_common,
                          CASSIA_parameters pCASSIA_parameters,
                          CASSIA_ratios pCASSIA_ratios,
                          Rcpp::DataFrame pCASSIA_sperling,
                          std::vector<double> extras_sperling,

                          double CH,
                          double B0,

                          bool LH_estim,
                          bool LN_estim,
                          bool LD_estim,

                          growth_values_out growth_in,
                          double last_year_HH,
                          int no_day);


/*
 *  Cell enlargement
 */

xylogenesis_out xylogenesis(int no_day,
                           int day,
                           CASSIA_parameters parameters,
                           CASSIA_common common,
                           bool xylogenesis_option,
                           bool environmental_effect_xylogenesis,
                           double TAir,

                           double n_rows,
                           double max_ew_cells,

                           double n_E_pot_old,
                           double n_W_pot_old,
                           double n_M_pot_old,

                           double g,
                           std::vector<double> en_growth_vector,

                           // Iterations
                           //double wall_tot,
                           //double tau_E,
                           //double tau_W,
                           //double n_E_old,
                           //double n_W_old,
                           //double n_M_old,
                           double tau_W_old,

                           // Parameters
                           double carbon_daily_rate_ew, // TODO: where does this number come from
                           double carbon_daily_rate_lw);

/*
 * Preles
 */

struct photosynthesis_out {
  double GPP;
  double ET;
  double SoilWater;
  double fS;
};

struct photo_out_vector {
  std::vector<double> GPP;
  std::vector<double> ET;
  std::vector<double> SoilWater;
  std::vector<double> fS;
};

photosynthesis_out preles_cpp(int day,
                              double I,
                              double T,
                              double P,
                              double D,
                              double CO2,
                              double fAPAR,
                              p1 Site_par,
                              p2 GPP_par,
                              p3 ET_par,
                              p4 SnowRain_par,
                              p5 Initials_snow,
                              double LOGFLAG);

Rcpp::DataFrame preles_test(Rcpp::DataFrame weather);

/*
 * Respiration
 */

respiration_out respiration(int day,
                            CASSIA_parameters parameters,
                            CASSIA_ratios ratios,
                            repola_out repola,
                            double TAir,
                            double TSoil,
                            bool temp_rise,
                            bool Rm_acclimation,
                            bool mN_varies,

                            // parameters that I am not sure about
                            double B0);

Rcpp::List respiration_test_cpp(Rcpp::DataFrame pCASSIA_parameters,
                                Rcpp::DataFrame pCASSIA_common,
                                Rcpp::DataFrame pCASSIA_ratios,
                                Rcpp::DataFrame pCASSIA_sperling,
                                std::vector<double> extras_sperling,
                                int ndays,
                                int day,
                                double TAir,
                                double TSoil,
                                bool temp_rise,
                                bool Rm_acclimation,
                                bool mN_varies,
                                double B0);

/*
 * Repola
 */

repola_out repola(CASSIA_parameters parameters);

/*
 * Actual growth
 */

growth_out actual_growth(CASSIA_parameters parameters,
                         CASSIA_common common,
                         carbo_tracker storage,
                         growth_out potential_growth,
                         respiration_out resp,
                         bool sperling_sugar_model);

/*
 * CASSIA_yearly
 */

Rcpp::List CASSIA_yearly(int start_year,
                         int end_year,

                         Rcpp::DataFrame weather,
                         std::vector<double> GPP_ref,

                         std::vector<double> pPREL,
                         Rcpp::DataFrame pCASSIA_parameters,
                         Rcpp::DataFrame pCASSIA_common,
                         Rcpp::DataFrame pCASSIA_ratios,
                         Rcpp::DataFrame pCASSIA_sperling,

                         double needle_mass_in, // The value of this should be 0 if you want the needle value to be calculated
                         double Throughfall,

                         Rcpp::List settings);

/*
 * Ecoevolutionary
 */

Rcpp::List CASSIA_eeo(int start_year,
                      int end_year,

                      Rcpp::DataFrame weather,
                      std::vector<double> GPP_ref,

                      std::vector<double> pPREL,
                      Rcpp::DataFrame pCASSIA_parameters,
                      Rcpp::DataFrame pCASSIA_common,
                      Rcpp::DataFrame pCASSIA_ratios,
                      Rcpp::DataFrame pCASSIA_sperling,
                      std::vector<double> parameters_R,

                      double needle_mass_in, // The value of this should be 0 if you want the needle value to be calculated
                      double Throughfall,
                      int trenching_year,

                      Rcpp::List settings);

/*
 * Soil
 */

Rcpp::List CASSIA_soil(int start_year,
                       int end_year,

                       Rcpp::DataFrame weather,
                       std::vector<double> GPP_ref,

                       std::vector<double> pPREL,
                       Rcpp::DataFrame pCASSIA_parameters,
                       Rcpp::DataFrame pCASSIA_common,
                       Rcpp::DataFrame pCASSIA_ratios,
                       Rcpp::DataFrame pCASSIA_sperling,
                       std::vector<double> parameters_R,

                       double needle_mass_in, // The value of this should be 0 if you want the needle value to be calculated
                       double Throughfall,
                       int trenching_year,

                       Rcpp::List settings);

/*
 * Parameters
 */

CASSIA_parameters make_CASSIA_parameters(Rcpp::DataFrame input_parameters, Rcpp::DataFrame input_sperling);
CASSIA_common make_common(Rcpp::DataFrame input);
CASSIA_ratios make_ratios(Rcpp::DataFrame input);
MYCOFON_function_out MYCOFON_structure_conversion(Rcpp::List input);


// Read from csv
CASSIA_ratios read_ratios(const std::string& filename, const std::string& site);
int CASSIA_ratios_test(std::string fratios, std::string site);
CASSIA_common read_common_parameters(const std::string& filename);
int CASSIA_common_test(std::string fcommon);
int CASSIA_parameter_test(CASSIA_parameters structure_params);

p1 make_p1(std::vector<double> input);
p2 make_p2(std::vector<double> input);
p3 make_p3(std::vector<double> input);
p4 make_p4(std::vector<double> input);
p5 make_p5(std::vector<double> input);
p7 make_p7(std::vector<double> input);

/*
 * Ring Width
 */

struct ring_width_out {
  double n_E_tot;
  double n_W_tot;
  double n_M_tot;

  double ew_cells_tot;

  double tot_mm;

  double max_ew_cells_tot;
  double max_ew_width_tot;
};

ring_width_out ring_width_generator(int day,
                                    ring_width_out previous_value,
                                    growth_values_out growth_previous,
                                    CASSIA_parameters parameters,
                                    double GD_tot);

/*
 * Settings defined
 */

// Define a struct to hold the settings
struct Settings {
  bool storage_reset;
  bool storage_grows;

  bool LN_estim;
  bool mN_varies;

  bool LD_estim;
  bool sD_estim_T_count;

  bool LH_estim;
  bool trees_grow;
  bool growth_decreases;
  bool needle_mass_grows;

  bool phloem_trigger;
  bool mycorrhiza;
  bool root_as_Ding;

  bool sperling_model;
  bool myco_model;
  bool xylogensis_option;

  bool PRELES_GPP;
  bool environmental_effect_xylogenesis;

  bool photosynthesis_as_input;
  bool phydro;
  bool fAPAR_Tian;

  int photoparameters;
  bool temp_rise;
  bool drought;
  bool Rm_acclimation;

  bool CASSIA_graphs;
  bool tests;
  bool etmodel;
  bool LOGFLAG;
};

// Inporting settings function
Settings parseSettings(Rcpp::List settingsList);

/*
 * Leap year
 */

int leap_year(int year);

/*
 * Plant Fate Logic
 */

//PlantAssimilationResult calc_plant_assimilation_rate(double fipar,
  //                                                   double PAR, double TAir, double VPD, double Precip, double CO2, double Nitrogen, double PA, double SWP,
    //                                                 phydro_canopy_parameters par, double lai, double n_layers, double crown_area, double height, double zeta);

#endif
