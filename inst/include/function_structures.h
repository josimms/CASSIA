#include <iostream>
#include <vector>

#ifndef FUNCTION_STRUCTURES_H
#define FUNCTION_STRUCTURES_H

/*
 * Respiration
 */

struct respiration_out {
  double RmN;
  double RmS;
  double RmR;
  double Rm_a;
};

struct resp_vector {
  std::vector<double> growth;
  std::vector<double> maintenance;
  std::vector<double> mycorrhiza;
  std::vector<double> microbes_FOM;
  std::vector<double> microbes_SOM;
};

/*
 * Repola
 */

struct repola_out {
  double needle_mass;
  double m_N;
  double m_N_tot;
  double m_R_tot;
};

/*
 * Growth
 */

struct growth_values_out {
  double sH;
  double fH;
  double HH;
  double sN;
  double fN;
  double sD;
  double fD;
  double sR;
  double fR;
  double n_rows;
  double GH;
  double GN;
  double GD;
  double max_N;
  double S_GPP;
  double dS_GPP;
  double S_GPP_ref;
  double dS_GPP_ref;
  double GPP_ref;
  double ew_cells_pot_max;
  double en_pot_growth;
  double pot_mm_max;
  double pot_mm;
  double wall_pot_growth;
  double n_E_pot;
  double n_W_pot;
  double n_M_pot;
  double n_E;
  double tau_E;
  double tau_W;
};

struct growth_out {
  double diameter;
  double needles;
  double height;
  double wall;
  double roots;
  double GD;
  double ecto;
  double bud;
  double use;
  double release;
  double respiration_growth;
  double respiration_maintenance;
  double g;
  growth_values_out previous_values;
};

struct growth_vector {
  std::vector<double> height;
  std::vector<double> needles;
  std::vector<double> roots;
  std::vector<double> diameter;
  std::vector<double> diameter_potential;
  std::vector<double> bud;
  std::vector<double> mycorrhiza;
  std::vector<double> ring_width;
  std::vector<double> height_tot;
  std::vector<double> wall_tot;
  std::vector<double> g;
  std::vector<double> en_pot_growth;
  std::vector<double> xylem_sh;
  std::vector<double> xylem_st;
  std::vector<double> phloem;
  std::vector<bool> tree_alive;
};

/*
 * Biomass
 */

struct biomass_out {
  double diameter;
  double needles;
  double height;
  double wall;
  double roots;
  double bud;
};

struct biomass_vector {
  std::vector<double> height;
  std::vector<double> needles;
  std::vector<double> roots;
  std::vector<double> diameter;
  std::vector<double> bud;
};

/*
 * Sugar model
 */

struct carbo_tracker_vector
{
  std::vector<double> needles; // Initial value from inputs
  std::vector<double> phloem; // Initial value from inputs
  std::vector<double> xylem_sh; // Initial value from inputs
  std::vector<double> xylem_st; // Initial value from inputs
  std::vector<double> roots; // Initial value from inputs
  std::vector<double> mycorrhiza; // Initial value 0
  double initial_amount; // Initial value from inputs
  double B;
};

struct conc_gradient
{
  double needles_to_phloem;
  double phloem_to_xylem_sh;
  double phloem_to_xylem_st;
  double phloem_to_roots;
  double roots_to_myco;
  double ratio_needles_to_phloem;
  double ratio_phloem_to_xylem_sh;
  double ratio_phloem_to_xylem_st;
  double ratio_phloem_to_roots;
};

typedef struct
{
  double needles; // Initial value from inputs
  double phloem; // Initial value from inputs
  double xylem_sh; // Initial value from inputs
  double xylem_st; // Initial value from inputs
  double roots; // Initial value from inputs
  double mycorrhiza; // Initial value 0
  double respiration;
  double surplus;
  double initial_amount; // Initial value from inputs
  double B;
} carbo_tracker;

struct carbo_values_out
{
  carbo_tracker storage_term;
  carbo_tracker Ad;
  carbo_tracker As;
  bool tree_alive;
  double sB0;
};

struct uptake_structre
{
  double ectomycorrhizal_transfer;
  double root_upatke;
  double ectomycorrhizal_upatke;
  double total_uptake;
};

struct uptake_structre_vector
{
  std::vector<double> ectomycorrhizal_transfer;
  std::vector<double> root_upatke;
  std::vector<double> ectomycorrhizal_upatke;
  std::vector<double> total_uptake;
};

struct carbo_balance
{
  carbo_tracker sugar;
  carbo_tracker starch;
  carbo_tracker storage;
  uptake_structre uptake;
  growth_out nitrogen_capacity;
  double resp_main;
  double resp_growth;
  carbo_values_out previous_values;
  double nitrogen_balance;
};

struct interaction
{
  double needles_to_phloem;
  double phloem_to_xylem_sh;
  double phloem_to_xylem_st;
  double phloem_to_roots;
};

struct sugar_values_vector
{
  std::vector<double> sugar;
  std::vector<double> starch;
  std::vector<double> storage;
  std::vector<double> starch_needles;
  std::vector<double> starch_phloem;
  std::vector<double> starch_xylem_sh;
  std::vector<double> starch_xylem_st;
  std::vector<double> starch_roots;
  std::vector<double> nitrogen_balance;
  std::vector<double> sugar_needles;
  std::vector<double> sugar_phloem;
  std::vector<double> sugar_xylem_sh;
  std::vector<double> sugar_xylem_st;
  std::vector<double> sugar_roots;
  std::vector<double> sugar_mycorrhiza;
  std::vector<double> sugar_surplus;
  std::vector<double> storage_needles;
  std::vector<double> storage_phloem;
  std::vector<double> storage_xylem_sh;
  std::vector<double> storage_xylem_st;
  std::vector<double> storage_roots;
  std::vector<double> storage_mycorrhiza;
  std::vector<double> nitrogen_capacity_needles;
  std::vector<double> nitrogen_capacity_bud;
  std::vector<double> nitrogen_capacity_wall;
  std::vector<double> nitrogen_capacity_height;
  std::vector<double> nitrogen_capacity_roots;
};

/*
 * Xylogenesis model
 */

struct xylo_previous_values {
  double sD;
  // double cell_d_ew; TODO: this is
  double cell_l;
  double M_suc;
  double gas_const;
  double DE_ew;
  double DE_lw;
  double carbon_daily_rate_ew;
  double carbon_daily_rate_lw;
  double n_rows;
  double tot_cells;
};

struct xylogensis_out {
  xylo_previous_values previous_value;
};

struct xylogenesis_out {
  double n_E;
  double n_E_pot;
  double n_W_pot;
  double n_M_pot;
};

/*
 * Needles
 */

struct needle_cohorts {
  double year_1;
  double year_2;
  double year_3;
};

/*
 * MYCOFON vector
 */

struct MYCOFON_vector {
  std::vector<double> C_biomass;
  std::vector<double> C_roots;
  std::vector<double> C_fungal;
  std::vector<double> N_roots;
  std::vector<double> N_fungal;
  std::vector<double> Respiration;
  std::vector<double> C_roots_NonStruct;
  std::vector<double> C_fungal_NonStruct;
  std::vector<double> N_roots_NonStruct;
  std::vector<double> N_fungal_NonStruct;
  std::vector<double> uptake_plant;
  std::vector<double> uptake_NH4_plant;
  std::vector<double> uptake_NO3_plant;
  std::vector<double> uptake_Norg_plant;
  std::vector<double> uptake_fungal;
  std::vector<double> uptake_NH4_fungal;
  std::vector<double> uptake_NO3_fungal;
  std::vector<double> uptake_Norg_fungal;
  std::vector<double> from_CASSIA;
  std::vector<double> to_CASSIA;
  std::vector<double> exudes_fungal;
  std::vector<double> exudes_plant;
  std::vector<double> Plant_demand;
  std::vector<double> Fungal_demand;
  std::vector<double> Plant_given;
  std::vector<double> Fungal_given;
};

/*
 * Gpp function output
 */

struct gpp_out {
  double fW;
  double fE;
  double fN;
  double fD;
  double fCO2;
  double gpp;
  double gpp380;
};

/*
 * phydro output
 */

struct PlantAssimilationResult{
  double gpp = 0;          ///< Gross plant-level production [kg-biomass yr-1]
  double npp = 0;          ///< Net plant-level production [kg-biomass yr-1]
  double trans = 0;        ///< Transpiration [kg-h2o yr-1]

  double dpsi_avg = 0;     ///< Soil-leaf water potential difference \f$\Delta\psi\f$
  double vcmax_avg = 0;    ///< Crown-area weighted average Vcmax across canopy layers [umol m-2 s-1]
  double vcmax25_avg = 0;  ///< Average Vcmax at 25 degC [umol m-2 s-1]
  double mc_avg = 0;       ///< \f$m_c = (\chi c_a - \Gamma^*)/(\chi c_a + K_M) \f$
  double gs_avg = 0;       ///< Crown-area weighted average stomatal conductance across canopy layers
  double c_open_avg = 0;   ///< Crown-area weighted average canopy opennness experience by the plant

  double nitrogen_avg = 0;       ///< Crown-area weighted average of Nitrogen across canopy laters

  double rleaf = 0;        ///< Leaf dark respiration rate [kg-biomass yr-1]
  double rroot = 0;        ///< Fine root respiration rate [kg-biomass yr-1]
  double rstem = 0;        ///< Sapwood respiration rate (excluding coarse root) [kg-biomass yr-1]

  double tleaf = 0;        ///< Leaf turnover rate [kg-biomass yr-1]
  double troot = 0;        ///< Fine root turnover rate [kg-biomass yr-1]
};

#endif
