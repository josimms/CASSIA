#include "CASSIA_types.h"
// CASSIA types includes function_structures.h, Rcpp.h, RcppCommon.h
// isostrem and vector are also in function_structures
#include "mycomodel.h"
// This includes soil and mycorrhizal parameters
#include "state_out.h"
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
#include <phydro.h>
#include <map>
#include <string>

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
};

struct eternal_parameters {
  double absolute_zero; // TOOD; what should this value be?
};

yearly_in yearly_initial_conditions(double days_per_year);


/*
 * sugar_model.cpp
 */

void sugar_model(int year,
                 int days_gone,
                 int day,
                 double TAir,
                 double PAR,
                 double PF,

                 const CASSIA_common& common,
                 const CASSIA_parameters& parameters,

                 double D00,
                 const growth_state& tree_state,

                 double& nitrogen_balance,
                 bool nitrogen_change,
                 bool nitrogen_contrast,

                 Settings boolsettings,

                 bool& tree_alive,
                 bool surplus_c,

                 double needles_mass, // Repola

                 carbo_tracker& sugar,
                 carbo_tracker& starch,
                 carbo_tracker& storage_term,
                 growth_out& nitrogen_capacity,
                 output_vector& out);

/*
 * growth.cpp
 */

void growth(int day,
            int days_gone,
            int year,

            growth_state& state,
            output_vector& out,

            double TAir,
            double TSoil_A,
            double TSoil_B,
            double Soil_Moisture,
            double PF,
            double GPP_ref,
            Settings boolsettings,

            const CASSIA_common& common,
            const CASSIA_parameters& parameters,
            const CASSIA_ratios& ratio,

            double CH,
            double B0,
            double GPP_mean,
            double GPP_previous_sum,

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
                              double LOGFLAG,
                              int CO2model);

Rcpp::DataFrame preles_test(Rcpp::DataFrame weather);

/*
 * Respiration
 */

void respiration(int day,
                 growth_state& tree_state,
                 output_vector& out,
                 const CASSIA_parameters& parameters,
                 const CASSIA_ratios& ratios,
                 double TAir,
                 double TSoil,
                 const Settings& settings,
                 double B0);

/*
 * Repola
 */

repola_out repola(CASSIA_parameters parameters);

/*
 * Actual growth
 */

void actual_growth(int day,
                   int days_gone,
                   const CASSIA_parameters& parameters,
                   const CASSIA_common& common,
                   const carbo_tracker& storage,
                   const photosynthesis_out& photosynthesis,
                   growth_state& tree_state,
                   output_vector& all_out,
                   Settings boolsettings,
                   growth_out nitrogen_capacity);

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

                         bool surplus_c,
                         bool nitrogen_change,
                         bool nitrogen_contrast,

                         double nitrogen_capacity,

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

phydro_canopy_parameters parPhydro_initalise(std::vector<double> phydro_params);

void print_phydro_parameters(const phydro_canopy_parameters& params);

/*
 * Ring Width
 */

void ring_width_generator(int day,
                          int days_gone,
                          growth_state& state,
                          output_vector& all_out,
                          const CASSIA_parameters& parameters);

// Inporting settings function
Settings parseSettings(Rcpp::List settingsList);

/*
 * Leap year
 */

int leap_year(int year);

/*
 * Plant Fate Logic
 */

PlantAssimilationResult calc_plant_assimilation_rate(double PAR, double PAR_max, double TAir, double VPD, double Precip, double CO2, double Nitrogen, double PA, double SWP,
                                                     phydro_canopy_parameters par, double lai, double crown_area, double height, double zeta, int day);

phydro::PHydroResultNitrogen leaf_assimilation_rate(double fipar, double fapar,
                                                    double PAR, double PAR_max, double TAir, double VPD, double Precip, double CO2, double Nitrogen, double PA, double SWP,
                                                    double TAir_assim, double PAR_assim, double VPD_assim, double CO2_assim, double SWP_assim, double PA_assim,
                                                    phydro_canopy_parameters par, double zeta);

void set_forcing_acclim(double TAir, double PAR, double VPD, double CO2, double SWP, double PA,
                        double& TAir_assim, double& PAR_assim, double& VPD_assim, double& CO2_assim, double& SWP_assim, double& PA_assim,
                        phydro_canopy_parameters par);

#endif


#ifndef READ_WEATHER_VARIABLES_H
#define READ_WEATHER_VARIABLES_H


struct weather_all {
  std::vector<double> TAir;
  std::vector<double> TSoil_A;
  std::vector<double> TSoil_B;
  std::vector<double> Soil_Moisture;
  std::vector<double> Precip;
  std::vector<double> Photosynthesis_IN;
  std::vector<double> Nitrogen;
  std::vector<double> PAR;
  std::vector<double> CO2;
  std::vector<double> VPD;
  std::vector<double> fAPAR;
  std::vector<double> PAR_max;
  std::vector<double> PA;
  std::vector<double> SWP;
};

weather_all readWeatherVariables(const Rcpp::DataFrame& weather, bool spp, bool preles, bool phydro);

#endif

/*
 * Initalise
 */

carbo_tracker carbo_tracker_init();

carbo_balance carbo_balance_init();

/*
 * Double to vectors
 */

void log_potential_growth(int day,
                          int days_gone,
                          const growth_state& tree_state,
                          output_vector& growth_out);

void log_actual_growth(int day,
                       int days_gone,
                       const growth_state& tree_state,
                       output_vector& growth_out);

void log_sugar(int day,
               int days_gone,
               const carbo_tracker& sugar,
               const carbo_tracker& starch,
               const carbo_tracker& storage_term,
               const growth_out& nitrogen_capacity,
               const double respiration_growth,
               const double respiration_maintenance,
               const double nitrogen_balance,
               const uptake_structre& uptake,
               const bool& tree_alive,
               output_vector& out);

void log_photosynthesis(int day,
                        int days_gone,
                        const photosynthesis_out& input,
                        output_vector& out);


/*
 * LAI and fAPAR
 */

void compute_fAPAR_used(int day,
                        int days_gone,
                        double LAI,
                        double max_needles,
                        const Settings& boolsettings,
                        const CASSIA_parameters& parameters,
                        const growth_state& tree_state,
                        const photosynthesis_out& photosynthesis,
                        const weather_all& climate,
                        output_vector& all_out);
/*
 * Make dataframes
 */

Rcpp::DataFrame createGrowthDataFrame(const output_vector& out);

// Generate DataFrame for sugar/starch/storage tracking
Rcpp::DataFrame createSugarDataFrame(const output_vector& out);

// Generate DataFrame for photosynthesis/preles outputs
Rcpp::DataFrame createPrelesDataFrame(const output_vector& out);

// Generate DataFrame for cumulative (culm) growth and related outputs
Rcpp::DataFrame createCulmGrowthDataFrame(const output_vector& out);

/*
 * Initalise
 */

void initialize_output_vector(output_vector& out, int simulation_time);

/*
 * Test
 */

void printOutputVectorSizes(const output_vector& out);
