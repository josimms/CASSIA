#include <iostream>
#include <vector>

#ifndef FUNCTION_STRUCTURES_H
#define FUNCTION_STRUCTURES_H

/*
 * Respiration
 */

struct repola_out {
  double needle_mass = 0.0;
  double m_N = 0.0;
  double m_N_tot = 0.0;
  double m_R_tot = 0.0;
};

struct resp_vector {
  std::vector<double> growth = {0.0};
  std::vector<double> maintenance = {0.0};
  std::vector<double> mycorrhiza = {0.0};
  std::vector<double> microbes_FOM = {0.0};
  std::vector<double> microbes_SOM = {0.0};

  void resize(size_t new_size) {
    growth.assign(new_size, 0.0);
    maintenance.assign(new_size, 0.0);
    mycorrhiza.assign(new_size, 0.0);
    microbes_FOM.assign(new_size, 0.0);
    microbes_SOM.assign(new_size, 0.0);
  }
};


/*
 * Growth state
 */

struct growth_state {
  // Ring width
  double n_E_tot = 0.0;
  double n_W_tot = 0.0;
  double n_M_tot = 0.0;
  double ew_cells_tot = 0.0;
  double lw_cells_tot = 0.0;
  double tot_mm = 0.0;
  double max_ew_cells_tot = 0.0;
  double max_ew_width_tot = 0.0;

  // Growth
  double sH = 0.0;
  double fH = 0.0;
  double HH = 0.0;
  double sN = 0.0;
  double fN = 0.0;
  double sD = 0.0;
  double fD = 0.0;
  double sR = 0.0;
  double fR = 0.0;
  double n_rows = 0.0;
  double GH = 0.0;
  double GN = 0.0;
  double GD = 0.0;
  double max_N = 0.0;
  double S_GPP = 0.0;
  double dS_GPP = 0.0;
  double S_GPP_ref = 0.0;
  double dS_GPP_ref = 0.0;
  double GPP_ref = 0.0;
  double ew_cells_pot_max = 0.0;
  double en_pot_growth = 0.0;
  double pot_mm_max = 0.0;
  double pot_mm = 0.0;
  double n_E_pot = 0.0;
  double n_W_pot = 0.0;
  double n_M_pot = 0.0;
  double n_E = 0.0;
  double tau_E = 0.0;
  double tau_W = 0.0;
  double wall_pot_growth = 0.0;

  // Growth out
  double g = 0.0;
  double diameter = 0.0;
  double needles = 0.0;
  double height = 0.0;
  double wall = 0.0;
  double roots = 0.0;
  double ecto = 0.0;
  double bud = 0.0;
  double use = 0.0;
  double release = 0.0;
  double ring_width = 0.0;

  // Respiration
  double RmN = 0.0;
  double RmS = 0.0;
  double RmR = 0.0;
  double Rm_a = 0.0;
};

/*
 * Growth
 */

struct growth_out {
  double diameter = 0.0;
  double needles = 0.0;
  double height = 0.0;
  double wall = 0.0;
  double roots = 0.0;
  double GD = 0.0;
  double ecto = 0.0;
  double bud = 0.0;
  double use = 0.0;
  double release = 0.0;
  double respiration_growth = 0.0;
  double respiration_maintenance = 0.0;
  double g = 0.0;
};

struct growth_vector {
  std::vector<double> height = {0.0};
  std::vector<double> needles = {0.0};
  std::vector<double> roots = {0.0};
  std::vector<double> diameter = {0.0};
  std::vector<double> diameter_potential = {0.0};
  std::vector<double> bud = {0.0};
  std::vector<double> mycorrhiza = {0.0};
  std::vector<double> ring_width = {0.0};
  std::vector<double> height_tot = {0.0};
  std::vector<double> wall_tot = {0.0};
  std::vector<double> g = {0.0};
  std::vector<double> en_pot_growth = {0.0};
  std::vector<double> xylem_sh = {0.0};
  std::vector<double> xylem_st = {0.0};
  std::vector<double> phloem = {0.0};
  std::vector<double> sapwood = {0.0};
  std::vector<double> leaf_mass = {0.0};
  std::vector<double> leaf_area = {0.0};
  std::vector<bool> tree_alive = {true};

  void resize(size_t new_size) {
    height.assign(new_size, 0.0);
    needles.assign(new_size, 0.0);
    roots.assign(new_size, 0.0);
    diameter.assign(new_size, 0.0);
    diameter_potential.assign(new_size, 0.0);
    bud.assign(new_size, 0.0);
    mycorrhiza.assign(new_size, 0.0);
    ring_width.assign(new_size, 0.0);
    height_tot.assign(new_size, 0.0);
    wall_tot.assign(new_size, 0.0);
    g.assign(new_size, 0.0);
    en_pot_growth.assign(new_size, 0.0);
    xylem_sh.assign(new_size, 0.0);
    xylem_st.assign(new_size, 0.0);
    phloem.assign(new_size, 0.0);
    sapwood.assign(new_size, 0.0);
    leaf_mass.assign(new_size, 0.0);
    leaf_area.assign(new_size, 0.0);
    tree_alive.assign(new_size, true);
  }

  void print(int day, const std::string& name = "growth_vector") const {
    std::cout << "\n=== " << name << " (Day " << day << ") ===\n";

    auto print_value = [&](const auto& vec, const std::string& label) {
      if (day >= 0 && static_cast<size_t>(day) < vec.size()) {
        std::cout << label << ": " << vec[day] << "\n";
      } else {
        std::cout << label << ": [Out of range index: " << day << "]\n";
      }
    };

    print_value(height, "height");
    print_value(needles, "needles");
    print_value(roots, "roots");
    print_value(diameter, "diameter");
    print_value(diameter_potential, "diameter_potential");
    print_value(bud, "bud");
    print_value(mycorrhiza, "mycorrhiza");
    print_value(ring_width, "ring_width");
    print_value(height_tot, "height_tot");
    print_value(wall_tot, "wall_tot");
    print_value(g, "g");
    print_value(en_pot_growth, "en_pot_growth");
    print_value(xylem_sh, "xylem_sh");
    print_value(xylem_st, "xylem_st");
    print_value(phloem, "phloem");
    print_value(sapwood, "sapwood");
    print_value(leaf_mass, "leaf_mass");
    print_value(leaf_area, "leaf_area");
    if (day >= 0 && static_cast<size_t>(day) < tree_alive.size()) {
      std::cout << "tree_alive: " << (tree_alive[day] ? "true" : "false") << "\n";
    } else {
      std::cout << "tree_alive: [Out of range index: " << day << "]\n";
    }

    std::cout << "=====================================\n";
  }
};


/*
 * Sugar model
 */

struct carbo_tracker_vector {
  std::vector<double> needles = {0.0}; // Initial value from inputs
  std::vector<double> phloem = {0.0}; // Initial value from inputs
  std::vector<double> xylem_sh = {0.0}; // Initial value from inputs
  std::vector<double> xylem_st = {0.0}; // Initial value from inputs
  std::vector<double> roots = {0.0}; // Initial value from inputs
  std::vector<double> to_mycorrhiza = {0.0}; // Initial value 0
  std::vector<double> surplus = {0.0}; // Initial value 0
  double initial_amount = 0.0; // Initial value from inputs
  double B = 0.0;

  void resize(size_t new_size) {
    needles.assign(new_size, 0.0);
    phloem.assign(new_size, 0.0);
    xylem_sh.assign(new_size, 0.0);
    xylem_st.assign(new_size, 0.0);
    roots.assign(new_size, 0.0);
    to_mycorrhiza.assign(new_size, 0.0);
    surplus.assign(new_size, 0.0);
    // Note: initial_amount and B remain unchanged
  }

  void print(int day, const std::string& name = "carbo_tracker_vector") const {
    std::cout << "\n=== " << name << " (Day " << day << ") ===\n";

    auto print_value = [&](const std::vector<double>& vec, const std::string& label) {
      if (day >= 0 && static_cast<size_t>(day) < vec.size()) {
        std::cout << label << ": " << vec[day] << "\n";
      } else {
        std::cout << label << ": [Out of range index: " << day << "]\n";
      }
    };

    print_value(needles,       "needles");
    print_value(phloem,        "phloem");
    print_value(xylem_sh,      "xylem_sh");
    print_value(xylem_st,      "xylem_st");
    print_value(roots,         "roots");
    print_value(to_mycorrhiza,"to_mycorrhiza");
    print_value(surplus,       "surplus");

    std::cout << "initial_amount: " << initial_amount << "\n";
    std::cout << "B: " << B << "\n";
    std::cout << "=====================================\n";
  }

};


struct conc_gradient {
  double needles_to_phloem = 0.0;
  double phloem_to_xylem_sh = 0.0;
  double phloem_to_xylem_st = 0.0;
  double phloem_to_roots = 0.0;
  double roots_to_myco = 0.0;
  double ratio_needles_to_phloem = 0.0;
  double ratio_phloem_to_xylem_sh = 0.0;
  double ratio_phloem_to_xylem_st = 0.0;
  double ratio_phloem_to_roots = 0.0;
};

struct carbo_tracker {
  double needles = 0.0; // Initial value from inputs
  double phloem = 0.0; // Initial value from inputs
  double xylem_sh = 0.0; // Initial value from inputs
  double xylem_st = 0.0; // Initial value from inputs
  double roots = 0.0; // Initial value from inputs
  double mycorrhiza = 0.0; // Initial value 0
  double respiration = 0.0;
  double surplus = 0.0;
  double initial_amount = 0.0; // Initial value from inputs
  double B = 0.0;
};

struct carbo_values_out {
  carbo_tracker storage_term = {};
  carbo_tracker Ad = {};
  carbo_tracker As = {};
  bool tree_alive = false;
  double sB0 = 0.0;
};

struct uptake_structre {
  double ectomycorrhizal_transfer = 0.0;
  double root_uptake = 0.0;
  double ectomycorrhizal_uptake = 0.0;
  double total_uptake = 0.0;
};

struct uptake_structre_vector {
  std::vector<double> ectomycorrhizal_transfer = {0.0};
  std::vector<double> root_uptake = {0.0};
  std::vector<double> ectomycorrhizal_uptake = {0.0};
  std::vector<double> total_uptake = {0.0};

  void resize(size_t new_size) {
    ectomycorrhizal_transfer.assign(new_size, 0.0);
    root_uptake.assign(new_size, 0.0);
    ectomycorrhizal_uptake.assign(new_size, 0.0);
    total_uptake.assign(new_size, 0.0);
  }
};


struct carbo_balance {
  carbo_tracker sugar = {};
  carbo_tracker starch = {};
  carbo_tracker storage = {};
  uptake_structre uptake = {};
  growth_out nitrogen_capacity = {};
  double resp_main = 0.0;
  double resp_growth = 0.0;
  carbo_values_out previous_values = {};
  double nitrogen_balance = 0.0;
};

struct interaction {
  double needles_to_phloem = 0.0;
  double phloem_to_xylem_sh = 0.0;
  double phloem_to_xylem_st = 0.0;
  double phloem_to_roots = 0.0;
};

/*
 * Xylogenesis model
 */

struct xylo_previous_values {
  double sD = 0.0;
  // double cell_d_ew; TODO: this is
  double cell_l = 0.0;
  double M_suc = 0.0;
  double gas_const = 0.0;
  double DE_ew = 0.0;
  double DE_lw = 0.0;
  double carbon_daily_rate_ew = 0.0;
  double carbon_daily_rate_lw = 0.0;
  double n_rows = 0.0;
  double tot_cells = 0.0;
};

struct xylogensis_out {
  xylo_previous_values previous_value = {};
};

struct xylogenesis_out {
  double n_E = 0.0;
  double n_E_pot = 0.0;
  double n_W_pot = 0.0;
  double n_M_pot = 0.0;
};

/*
 * Needles
 */

struct needle_cohorts {
  double year_1 = 0.0;
  double year_2 = 0.0;
  double year_3 = 0.0;
};

/*
 * MYCOFON vector
 */

struct MYCOFON_vector {
  std::vector<double> C_biomass = {0.0};
  std::vector<double> C_roots = {0.0};
  std::vector<double> C_fungal = {0.0};
  std::vector<double> N_roots = {0.0};
  std::vector<double> N_fungal = {0.0};
  std::vector<double> Respiration = {0.0};
  std::vector<double> C_roots_NonStruct = {0.0};
  std::vector<double> C_fungal_NonStruct = {0.0};
  std::vector<double> N_roots_NonStruct = {0.0};
  std::vector<double> N_fungal_NonStruct = {0.0};
  std::vector<double> uptake_plant = {0.0};
  std::vector<double> uptake_NH4_plant = {0.0};
  std::vector<double> uptake_NO3_plant = {0.0};
  std::vector<double> uptake_Norg_plant = {0.0};
  std::vector<double> uptake_fungal = {0.0};
  std::vector<double> uptake_NH4_fungal = {0.0};
  std::vector<double> uptake_NO3_fungal = {0.0};
  std::vector<double> uptake_Norg_fungal = {0.0};
  std::vector<double> from_CASSIA = {0.0};
  std::vector<double> to_CASSIA = {0.0};
  std::vector<double> exudes_fungal = {0.0};
  std::vector<double> exudes_plant = {0.0};
  std::vector<double> Plant_demand = {0.0};
  std::vector<double> Fungal_demand = {0.0};
  std::vector<double> Plant_given = {0.0};
  std::vector<double> Fungal_given = {0.0};
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

/*
 * Settings defined
 */

// Define a struct to hold the settings
struct Settings {
  bool storage_reset = false;
  bool storage_grows = false;

  bool LN_estim = false;
  bool mN_varies = false;

  bool LD_estim = false;
  bool sD_estim_T_count = false;

  bool LH_estim = false;
  bool trees_grow = false;
  bool growth_decreases = false;
  bool needle_mass_grows = false;

  bool phloem_trigger = false;
  bool mycorrhiza = false;
  bool root_as_Ding = false;

  bool sperling_model = false;
  bool myco_model = false;
  bool xylogensis_option = false;

  bool environmental_effect_xylogenesis = false;

  bool photosynthesis_as_input = false;
  bool preles = false;
  bool phydro = false;
  bool fAPAR_Tian = false;

  int photoparameters = 0;
  bool temp_rise = false;
  bool drought = false;
  bool Rm_acclimation = false;

  bool CASSIA_graphs = false;
  bool tests = false;
  bool etmodel = false;
  bool LOGFLAG = false;
};

/*
 * Preles
 */

struct photosynthesis_out {
  double GPP = 0.0;
  double ET = 0.0;
  double SoilWater = 0.0;
  double fS = 0.0;
};

struct photo_out_vector {
  std::vector<double> GPP = {0.0};
  std::vector<double> ET = {0.0};
  std::vector<double> SoilWater = {0.0};
  std::vector<double> fS = {0.0};
  std::vector<double> fAPAR = {0.0};

  void resize(size_t new_size) {
    GPP.assign(new_size, 0.0);
    ET.assign(new_size, 0.0);
    SoilWater.assign(new_size, 0.0);
    fS.assign(new_size, 0.0);
    fAPAR.assign(new_size, 0.0);
  }
};

/*
 * Output vector
 */

struct output_vector {
  /*
   * Years and Days
   */
  std::vector<int> year = {0};
  std::vector<int> day = {0};

  /*
   * Growth
   */
  // Potential
  std::vector<double> potential_diameter = {0.0};
  std::vector<double> potential_needles = {0.0};
  std::vector<double> potential_height = {0.0};
  std::vector<double> potential_wall = {0.0};
  std::vector<double> potential_roots = {0.0};
  std::vector<double> potential_bud = {0.0};
  std::vector<double> potential_ecto = {0.0};
  std::vector<double> potential_ring_width = {0.0};

  // Actual
  std::vector<double> g = {0.0};
  std::vector<double> diameter = {0.0};
  std::vector<double> needles = {0.0};
  std::vector<double> height = {0.0};
  std::vector<double> wall = {0.0};
  std::vector<double> roots = {0.0};
  std::vector<double> ecto = {0.0};
  std::vector<double> bud = {0.0};
  std::vector<double> use = {0.0};
  std::vector<double> release = {0.0};

  std::vector<double> ring_width = {0.0};
  std::vector<double> en_pot_growth = {0.0};

  /*
   * Culmative Growth
   */
  growth_vector culm_growth = {};

  /*
   * Respiration
   */
  std::vector<double> RmN = {0.0};
  std::vector<double> RmS = {0.0};
  std::vector<double> RmR = {0.0};
  std::vector<double> Rm_a = {0.0};

  /*
   * Repola
   */
  std::vector<double> needle_mass = {0.0};
  std::vector<double> m_N = {0.0};
  std::vector<double> m_N_tot = {0.0};
  std::vector<double> m_R_tot = {0.0};

  std::vector<double> fAPAR = {0.0};
  std::vector<double> LAI = {0.0};

  /*
   * Sugar
   */
  carbo_tracker_vector sugar_vector = {};
  carbo_tracker_vector starch_vector = {};
  carbo_tracker_vector storage_term_vector = {};
  resp_vector respiration_output = {};
  growth_vector nitrogen_capacity_vector = {};
  std::vector<double> nitrogen_balance = {0.0};
  uptake_structre_vector uptake_vector = {};
  std::vector<bool> tree_alive = {true};
  std::vector<double> sugar = {0.0};
  std::vector<double> starch = {0.0};

  /*
   * Photosynthesis
   */
  photo_out_vector photosynthesis = {};
};

#endif
