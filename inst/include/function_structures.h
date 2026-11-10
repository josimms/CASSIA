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
  // Nitrogen retranslocated from today's needle/root/mycorrhiza litter
  // (actual_growth.cpp), read by sugar_model.cpp's nitrogen_balance
  // update on the NEXT day's call - the same one-day lag
  // uptake.total_uptake already relies on (sugar_model.cpp reads it at
  // the top of a call, before this same day's actual_growth() has run
  // yet to set tomorrow's value).
  double N_retranslocated_today = 0.0;

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

  // Makela & Valentine (2020) pipe-model correction, roots referenced
  // against LEAF MASS via Eq. 4.13 (W_r = alpha_r * W_f) - supersedes an
  // earlier leaf-area version (see development_notebooks/crown_height.Rmd
  // Sect. 5/6 and manuscript/paper_plan.Rmd's "Root growth" section for
  // why: Repola's own root equations, Eq. 10/11, predict root mass from
  // diameter alone, with no leaf link at all, so this pipe-model ratio is
  // the only formula in CASSIA that actually ties root mass to leaf mass).
  // Recomputed once per year (day == 0), not daily like height's, and
  // held here across the rest of that year's days.
  double C_R_year = 1.0;

  // Running max of roots across the CURRENT year, updated daily in
  // actual_growth.cpp; read and reset at day == 0 in growth.cpp, at
  // which point it holds the just-finished year's max (not yet reset for
  // the new year). C_R_year is computed from this and
  // max_leaf_mass_running_for_roots below, not the previous year's
  // last-day value - a within-year peak, not wherever the pool happens
  // to sit at the literal year-end snapshot (which can undershoot the
  // peak once senescence/litter are active).
  double max_roots_running = 0.0;

  // A SECOND running-max tracker for leaf_mass, dedicated to the root
  // correction above - not a simplification choice, a necessity: this
  // file's needle section (leaf_mass_correction's own day == 0 block,
  // growth.cpp) already resets max_leaf_mass_running (below) to 0 for the
  // new year BEFORE this file's root section runs, within the same
  // growth() call - reading that field from the root section would
  // always see 0, never last year's peak. Tracked unconditionally in
  // actual_growth.cpp (like max_roots_running), independent of
  // leaf_mass_correction, so it stays ready if root_leaf_correction is
  // toggled on alone - though in that case leaf_mass itself is the old
  // memoryless, annually-reset quantity (fAPAR_and_LAI.cpp), not a
  // genuinely accumulating pool, so this correction is only really
  // meaningful paired with leaf_mass_correction = TRUE.
  double max_leaf_mass_running_for_roots = 0.0;

  // This year's root balance target, alpha_r * (last year's leaf_mass
  // peak) - the same quantity growth.cpp's local R_bal variable computes,
  // threaded through to actual_growth.cpp so root litter's own boost
  // factor (makela_loss_boost_factor) can react to the same fixed,
  // externally-anchored target as growth's C_R_year, not the moving
  // roots pool itself. Mirrors N_repo_current's role for needles.
  double R_bal_current = 0.0;

  // Root/mycorrhiza turnover: same sticky-senescence mechanism as
  // fAPAR_and_LAI.cpp's fS_reached_one for needles, kept separately here
  // since that flag isn't itself passed through to actual_growth(). Reset
  // to false at day == 0.
  bool fS_reached_one_roots = false;

  // leaf_mass_correction (Makela correction on needle leaf_mass,
  // 2026-11-05 redesign): C_N_year damps state.needles DIRECTLY
  // (growth.cpp) - the real needle-growth demand signal sugar_model.cpp
  // already uses - rather than driving a separate parallel growth term.
  // Computed fresh each day against max_leaf_mass_running (THIS year's
  // own running PEAK so far, reset to 0 at day == 0, updated daily in
  // actual_growth.cpp): the peak's real-time rise-with-growth gives this
  // correction within-season negative feedback; a last-year-frozen or
  // winter-stock (running minimum) reference were both tried and found
  // worse (unstable / overshooting respectively - see growth.cpp).
  double C_N_year = 1.0;
  double max_leaf_mass_running = 0.0;

  // One-way latch (growth.cpp, reset false at day==0): true once this
  // year's needle growth (fN) has definitively ended (sN has crossed
  // sNc^2), not just temporarily 0 before the season starts. Used in
  // actual_growth.cpp to freeze max_leaf_mass_running once growth is
  // genuinely over, so litter's correction reacts to a settled peak
  // rather than one that could still be rising from residual overlap
  // between the growth season and litter's own day>=212 gate (confirmed
  // real in ~1/4 of simulated years, see growth.cpp's comment).
  bool needle_growth_season_ended = false;

  // Needle litter timing: a skew-normal density over day-of-year
  // (parameters.omega_litter/alpha_litter), replacing the earlier
  // fS^2-accumulator-plus-Makela-bell design (2026-11-10 notebook
  // comparison: validation RMSE 0.068 vs 0.153 for the old gated
  // mechanism - development_notebooks/leaf_mass_correction.Rmd Sect. 5.1).
  // sum_degree_fL_running/last_year are the same one-year-lagged
  // self-normalizer the old mechanism already used (freeze-at-day-0,
  // "only overwrite if running>epsilon" bootstrap-year guard) - this
  // part of the design is unchanged, only what gets accumulated differs.
  //
  // xi (the skew-normal's location) is itself one-year-lagged for the
  // same reason: the day fS first drops below 0.5 each year can only be
  // known once that year is already most of the way over, so THIS
  // year's xi is set from LAST year's crossing day (+ a fixed
  // parameters.xi_offset_litter), not this year's own.
  double sum_degree_fL_running = 0.0;
  double sum_degree_fL_last_year = 0.0;
  double fS_cross_day_this_year = -1.0;  // 1-indexed day-of-year, -1 = not yet found this year
  bool   fS_crossed_this_year = false;
  double fS_cross_day_last_year = -1.0;  // frozen at day==0; -1 = no prior year yet (bootstrap)
  double xi_litter = -1.0;               // this year's skew-normal location; -1 = invalid/bootstrap year

  // Needle age cohorts (2026-11-10), replacing the old CASSIA.cpp-local,
  // never-actually-used needle_cohorts struct ("NOTE: not used anywhere!").
  // One entry per age-class, sized to parameters.n_age at first use (3 for
  // Hyde, 5 for Vario) - [0] is the newest (this year's real growth,
  // frozen at the NEXT year's day==0), [n_age-1] is the oldest, about to
  // be shed. Empty until the first day==0 call seeds it (growth.cpp).
  std::vector<double> needle_cohort_mass;
  double needle_litter_magnitude_this_year = 0.0;  // oldest cohort's mass, read at day==0 before rollover - this year's needle_litter total

  // Cone litter (2026-11-10): a masting cycle, independent of the needle
  // timing above - cones take multiple years to mature in Scots pine, and
  // masting is highly irregular/resource-driven rather than a fitted
  // continuous shape (confirmed in leaf_mass_correction.Rmd Sect. 5.2:
  // RMSE-fitting a cone window just inflates duration chasing masting
  // noise). cone_biomass_standing accumulates every year from real bud
  // growth (out.culm_growth.bud, the same per-year running total pattern
  // needle_cohort_mass reads from out.culm_growth.needles); every
  // n_age-th year (years_since_masting resets to 0 - growth.cpp reuses
  // n_age rather than a second hardcoded cycle length, so the period is
  // site-dynamic: 3 Hyde, 5 Vario), the WHOLE accumulated stock sheds as
  // litter (cone_litter_magnitude_this_year) and standing drops back to
  // 0 - a mature cone is a discrete unit, unlike needles' continuous
  // age-class turnover, so nothing is retained fractionally.
  double cone_biomass_standing = 0.0;
  int    years_since_masting = 0;
  double cone_litter_magnitude_this_year = 0.0;  // 0 except masting years

  // This year's N_repo (growth.cpp), threaded through to actual_growth.cpp
  // for the litter calculation - litter is scaled against this fixed,
  // externally-anchored target (N_repo/n_age), not the moving leaf_mass
  // pool itself, so a pool that's already drifted above target doesn't
  // shed proportionally more just because it's bigger.
  double N_repo_current = 0.0;

  // Litter's own correction factor (makela_loss_boost_factor, growth.cpp)
  // is NOT stored here - unlike C_N_year, it's computed fresh each day in
  // actual_growth.cpp from max_leaf_mass_running directly, to avoid a
  // one-year-lagged feedback signal that was confirmed unstable (see
  // growth.cpp's comment).

  // Respiration
  double RmN = 0.0;
  double RmS = 0.0;
  double RmR = 0.0;
  double Rm_a = 0.0;

  // Shared seasonality-model state (Mäkelä et al. 2004's fS, Rprebasso's
  // fS_model()/gpp.c) - 2026-11-10. Previously only ever computed inside
  // preles_cpp()'s own static local, so leaf_mass_correction's fS^2
  // needle-litter timing and the root-senescence gate above
  // (fS_reached_one_roots) only worked under preles=TRUE; now computed
  // once per day in CASSIA.cpp's year loop regardless of which
  // photosynthesis model is active (see
  // fS_decoupled_from_preles_implementation_plan.md), and preles_cpp()
  // consumes this same value instead of keeping its own separate copy.
  // 20.0 matches PRELES's own Initials_snow.S baseline (preles_test.cpp)
  // - reset to this on the same condition preles_cpp()'s reset_state
  // already uses (day==0 && (soil_water_reset || year==start_year)), not
  // left continuous across years, so each spring's ramp-up starts fresh
  // from the same baseline PRELES always has.
  double fS_temp_sum = 20.0;
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
  std::vector<double> root_litter = {0.0};
  std::vector<double> mycorrhiza_litter = {0.0};
  std::vector<double> needle_litter = {0.0};
  std::vector<double> cone_litter = {0.0};
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
    root_litter.assign(new_size, 0.0);
    mycorrhiza_litter.assign(new_size, 0.0);
    needle_litter.assign(new_size, 0.0);
    cone_litter.assign(new_size, 0.0);
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
  // PRELES soil-water/snow/canopy/temperature-sum state (preles_test.cpp's
  // preles_cpp()) resets to Initials_snow's fixed values every year by
  // default (TRUE), matching storage_reset's own annual-reset semantics.
  // FALSE carries that state over continuously across years instead -
  // needed for ORCHIMIC's continuous multi-year runs (see
  // manscript_orchemic/PROJECT_INSTRUCTIONS.md). Added 2026-11-09.
  bool soil_water_reset = true;
  bool storage_grows = false;

  bool LN_estim = false;
  bool mN_varies = false;

  bool LD_estim = false;
  bool sD_estim_T_count = false;

  bool LH_estim = false;
  bool trees_grow = false;
  bool growth_decreases = false;
  bool needle_mass_grows = false;
  bool height_diameter_correction = false;
  // tau_H moved to CASSIA_parameters (2026-11-09) - see final_parameters.h.

  // Makela & Valentine (2020) pipe-model correction for root growth,
  // referenced against leaf mass via the pipe model's own ratio
  // (Ch.4, Eq. 4.13: W_r = alpha_r * W_f) - see growth.cpp's root
  // section and manuscript/paper_plan.Rmd's "Root growth" section.
  // Supersedes an earlier leaf-area-referenced version, 2026-11-06.
  bool root_leaf_correction = false;
  // tau_R moved to CASSIA_parameters (2026-11-09) - see final_parameters.h.

  // Root/mycorrhiza litter (turnover). When root_leaf_correction is also
  // TRUE, root_litter is additionally boosted/damped toward the same
  // Eq. 4.13 target (actual_growth.cpp); mycorrhiza_litter is not - no
  // pipe-model analogue exists for mycorrhiza.
  bool root_turnover = false;

  // Makela (1999) structural-balance correction for needle leaf_mass,
  // referenced against Repola's own needle-mass prediction (growth.cpp,
  // actual_growth.cpp) - replaces fAPAR_and_LAI.cpp's memoryless,
  // annually-reset leaf_mass formula with a genuinely accumulating pool
  // for THIS quantity only; leaf_area/sapwood/fAPAR stay on the old
  // pathway unchanged. See manuscript/paper_plan.Rmd Working Notes.
  bool leaf_mass_correction = false;
  // tau_N, sLc_litter, and theta_N_*/retrans_* (nitrogen content of
  // needle/root/mycorrhiza litter net of retranslocation) all moved to
  // CASSIA_parameters (2026-11-09), along with tau_H/tau_R above - see
  // final_parameters.h for the fields and data-raw/building_data.R for
  // sources/references. CASSIA_cpp()'s only non-toggle arguments are now
  // the big parameter/common/etc. tables.

  // Cone litter (masting cycle, growth.cpp/actual_growth.cpp) - its own
  // toggle (2026-11-10), split out from leaf_mass_correction. It was
  // originally added inside that mechanism's own code blocks purely
  // because that's where the needle work was happening, but it has no
  // real dependency on leaf_mass_correction's needle cohorts/skew-normal
  // - only on out.culm_growth.bud (always tracked, unconditional) and
  // its own independent state (cone_biomass_standing/years_since_masting/
  // cone_litter_magnitude_this_year). Same false default/opt-in pattern
  // as root_turnover above.
  bool cone_turnover = false;

  bool phloem_trigger = false;
  bool mycorrhiza = false;
  bool root_as_Ding = false;

  bool organ_level_sugar = false;
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

  bool soil_moisture_effect_on_shoot = false;
  bool soil_moisture_effect_on_needles = false;
  bool soil_moisture_effect_on_diameter = false;
  std::vector<int> driver_N = {1, 1};
  std::vector<int> driver_H = {0, 2};
  std::vector<int> driver_D = {1, 1};
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

  // Daily (not cumulative) litter flux - contrast with culm_growth's own
  // root_litter/mycorrhiza_litter, which are running totals. root_litter/
  // mycorrhiza_litter are 0 whenever root_turnover = FALSE; needle_litter
  // is unconditional (read directly off leaf_mass's already-existing
  // senescence decline, not behind a toggle - see fAPAR_and_LAI.cpp).
  std::vector<double> root_litter = {0.0};
  std::vector<double> mycorrhiza_litter = {0.0};
  std::vector<double> needle_litter = {0.0};

  // Cone litter (2026-11-10): masting-driven, see growth_state's
  // cone_biomass_standing/years_since_masting (function_structures.h
  // above) and actual_growth.cpp for the fixed day 121-273 shedding
  // window. 0 except masting years (every 3rd year). No dedicated N
  // model yet - cone_litter_N is a placeholder, always 0, same reasoning
  // as branch_litter_N below.
  std::vector<double> cone_litter = {0.0};
  std::vector<double> cone_litter_N = {0.0};

  // Nitrogen content of the ABOVE litter fluxes, net of retranslocation -
  // i.e. the N actually leaving the tree in shed tissue, not the gross N
  // it contained before resorption (actual_growth.cpp). Zero whenever the
  // corresponding biomass litter above is zero (root_turnover = FALSE,
  // etc.).
  std::vector<double> needle_litter_N = {0.0};
  std::vector<double> root_litter_N = {0.0};
  std::vector<double> mycorrhiza_litter_N = {0.0};

  // Always 0 - CASSIA has no standalone branch biomass pool yet (branches
  // are folded into xylem/sapwood form-factor accounting only, see
  // fAPAR_and_LAI.cpp). Placeholder output (mass and N) so a future
  // branch pool has somewhere to log litter without changing the
  // Turnover data frame's shape again.
  std::vector<double> branch_litter = {0.0};
  std::vector<double> branch_litter_N = {0.0};

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
