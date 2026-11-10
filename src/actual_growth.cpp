#include "CASSIA.h"

void actual_growth(int day,
                   int days_gone,
                   const CASSIA_parameters& parameters,
                   const CASSIA_common& common,
                   const carbo_tracker& storage,
                   const photosynthesis_out& photosynthesis,
                   growth_state& tree_state,
                   output_vector& all_out,
                   Settings boolsettings,
                   growth_out nitrogen_capacity) {

  // Nitrogen content of litter, net of retranslocation (2026-11-05).
  // N:C thresholds and retranslocation fractions are real, sourced
  // CASSIA_parameters fields (parameters.theta_N_*/retrans_*, see
  // data-raw/building_data.R for sourcing) rather than hardcoded here.
  double needle_litter_N_retrans = 0.0, root_litter_N_retrans = 0.0, myco_litter_N_retrans = 0.0;

  // CASSIA has no standalone branch biomass pool yet (see fAPAR_and_LAI.cpp's
  // form-factor-only branch split) - mass and N populated with 0 so the
  // Turnover data frame already has a slot for both once a pool exists.
  all_out.branch_litter[days_gone + day]   = 0.0;
  all_out.branch_litter_N[days_gone + day] = 0.0;

  double denom = all_out.culm_growth.phloem[day + days_gone] + all_out.culm_growth.xylem_sh[day + days_gone] + all_out.culm_growth.xylem_st[day + days_gone];

  double phloem_share   = denom > 0 ? all_out.culm_growth.phloem[day + days_gone] / denom : 0.0;
  double xylem_st_share = denom > 0 ? all_out.culm_growth.xylem_st[day + days_gone] / denom : 0.0;
  double xylem_sh_share = denom > 0 ? all_out.culm_growth.xylem_sh[day + days_gone] / denom : 0.0;
  /*
   * Height
   */

  double storage_height = 0;
  if (boolsettings.organ_level_sugar) {
    storage_height = (phloem_share * storage.phloem + xylem_st_share * storage.xylem_st + xylem_sh_share * storage.xylem_sh);
  } else {
    storage_height = storage.needles;
  }
  tree_state.height = tree_state.height * std::min(storage_height, nitrogen_capacity.height);

  /*
   * Wall
   */

  double storage_wall = 0;
  if (boolsettings.organ_level_sugar) {
    storage_wall = (phloem_share * storage.phloem + xylem_st_share * storage.xylem_st + xylem_sh_share * storage.xylem_sh);
  } else {
    storage_wall = storage.needles;
  }
  tree_state.wall = tree_state.diameter * std::min(storage_wall, nitrogen_capacity.wall);

  /*
   * Bud
   */
  tree_state.bud = tree_state.bud * std::min(storage.needles, nitrogen_capacity.bud);

  /*
   * Needles
   */
  tree_state.needles = tree_state.needles * std::min(storage.needles, nitrogen_capacity.needles);

  /*
   * Roots
   */
  double storage_roots;
  if (boolsettings.organ_level_sugar) {
    storage_roots = storage.roots;
  } else {
    storage_roots = storage.needles;
  }
  tree_state.roots = tree_state.roots * std::min(storage_roots, nitrogen_capacity.roots);

  /*
   * Mycorrhiza
   *
   * No dedicated nitrogen-capacity threshold exists for mycorrhiza (see
   * nitrogen_storage() in sugar_model.cpp); it is gated the same way as
   * roots, immediately above.
   */
  double storage_mycorrhiza;
  if (boolsettings.organ_level_sugar) {
    storage_mycorrhiza = storage.roots;
  } else {
    storage_mycorrhiza = storage.needles;
  }
  tree_state.ecto = tree_state.ecto * std::min(storage_mycorrhiza, nitrogen_capacity.roots);

  /*
   * GD
   */
  double storage_GD;
  if (boolsettings.organ_level_sugar) {
    storage_GD = (phloem_share * storage.phloem + xylem_st_share * storage.xylem_st + xylem_sh_share * storage.xylem_sh);
  } else {
    storage_GD = storage.needles;
  }
  tree_state.GD = tree_state.GD * std::min(storage_GD, nitrogen_capacity.wall);

  /*
   * Save the output to a vector
   */
  log_actual_growth(day, days_gone, tree_state, all_out);

  /*
   * Culmative Log Manual
   */
  int index_ref = days_gone + day - 1;
  if (index_ref < 0) {
    index_ref = 0;
  }

  // Culmulative
  all_out.culm_growth.height[days_gone + day]   = all_out.culm_growth.height[index_ref] + tree_state.height;

  // Root/mycorrhiza turnover (litter): same senescence SIGNAL as leaf
  // area's cohort turnover (fAPAR_and_LAI.cpp, Tian 2021) - photosynthesis.fS,
  // sticky once it first reaches 1.0 after day 50 (fS_reached_one_roots),
  // tracking fS's own decline from then on. Deliberately NOT a literal
  // copy of the needle formula though: that formula's flat (n_age-1)/n_age
  // daily decay only stays sensible because it is paired with an equally
  // large, phenology-driven regrowth term (f_modifer) that has no root
  // equivalent here - transplanting it as-is would drain the pool during
  // every zero-growth winter day. Using today's actual realized growth
  // (tree_state.roots/.ecto, already correct - see the previous commit)
  // as the replacement term instead, and turning the senescence signal
  // into a DEGREE of senescence (1 - fS) so loss is ~0 except during the
  // real autumn decline, scaled by a per-day rate derived from a
  // turnover-time placeholder (see root_turnover_years below).
  //
  // NOTE: parameters.root_lifetime looks like the natural rate to reuse
  // here (same role as n_age plays for needles) but is NA for this site's
  // parameter table - it is currently only read by growth.cpp's
  // root_as_Ding == FALSE branch, which nothing in this session's runs
  // exercises (root_as_Ding defaults TRUE). Using it as-is would silently
  // NaN the whole run the moment root_turnover is enabled, so a separate,
  // explicitly-flagged placeholder is used instead until root_lifetime
  // itself is actually populated/validated.
  if (day == 0) {
    tree_state.fS_reached_one_roots = false;
  }
  if (photosynthesis.fS >= 1.0 && day > 50) {
    tree_state.fS_reached_one_roots = true;
  }
  double root_senescence_degree = tree_state.fS_reached_one_roots ? (1.0 - photosynthesis.fS) : 0.0;

  if (boolsettings.root_turnover) {
    const double root_turnover_years = 3.0; // TODO: placeholder (matches n_age's default) - needs real calibration, and ideally parameters.root_lifetime once that's actually populated
    double daily_turnover_fraction = root_senescence_degree / (root_turnover_years * 365.0);

    // Pipe-model litter boost (mirrors needle litter's litter_boost_today,
    // growth.cpp's makela_loss_boost_factor): symmetric, uncapped ratio of
    // this year's running root peak (tree_state.max_roots_running) to this
    // year's fixed alpha_r*leaf_mass target (tree_state.R_bal_current, set
    // in growth.cpp's root_leaf_correction block) - > 1 sheds roots faster
    // once the pool has overgrown past the pipe-model target, < 1 sheds
    // slower while undergrown, same reasoning as needles' own correction.
    // Naturally stays 1.0 whenever root_leaf_correction has never run
    // (R_bal_current still 0 - makela_loss_boost_factor's own epsilon
    // guard), so this is safe to leave unconditional here. Scoped to roots
    // only, not mycorrhiza - Eq. 4.13 is a fine-root relationship (Makela
    // & Valentine 2020, Ch.4); mycorrhiza has no pipe-model analogue, so
    // myco_litter_today stays on the base senescence-driven fraction
    // alone, same as before this change.
    double root_bal_boost_today = makela_loss_boost_factor(
      tree_state.R_bal_current, tree_state.max_roots_running, parameters.tau_R);

    double root_litter_today = all_out.culm_growth.roots[index_ref] * daily_turnover_fraction * root_bal_boost_today;
    double myco_litter_today = all_out.culm_growth.mycorrhiza[index_ref] * daily_turnover_fraction;

    all_out.culm_growth.roots[days_gone+day]      = std::max(0.0, all_out.culm_growth.roots[index_ref]      - root_litter_today + tree_state.roots);
    all_out.culm_growth.mycorrhiza[days_gone+day] = std::max(0.0, all_out.culm_growth.mycorrhiza[index_ref] - myco_litter_today + tree_state.ecto);

    all_out.culm_growth.root_litter[days_gone+day]       = all_out.culm_growth.root_litter[index_ref]       + root_litter_today;
    all_out.culm_growth.mycorrhiza_litter[days_gone+day] = all_out.culm_growth.mycorrhiza_litter[index_ref] + myco_litter_today;

    // Daily (non-cumulative) flux - see "Turnover" output dataframe.
    all_out.root_litter[days_gone+day]       = root_litter_today;
    all_out.mycorrhiza_litter[days_gone+day] = myco_litter_today;

    double root_litter_N_gross = parameters.theta_N_roots / parameters.carbon_share * root_litter_today;
    double myco_litter_N_gross = parameters.theta_N_myco  / parameters.carbon_share * myco_litter_today;
    root_litter_N_retrans = root_litter_N_gross * parameters.retrans_roots;
    myco_litter_N_retrans = myco_litter_N_gross * parameters.retrans_myco;
    all_out.root_litter_N[days_gone+day]       = root_litter_N_gross - root_litter_N_retrans;
    all_out.mycorrhiza_litter_N[days_gone+day] = myco_litter_N_gross - myco_litter_N_retrans;
  } else {
    all_out.culm_growth.roots[days_gone+day]      = all_out.culm_growth.roots[index_ref]      + tree_state.roots;
    all_out.culm_growth.mycorrhiza[days_gone+day] = all_out.culm_growth.mycorrhiza[index_ref] + tree_state.ecto;

    all_out.root_litter[days_gone+day]       = 0.0;
    all_out.mycorrhiza_litter[days_gone+day] = 0.0;
    all_out.root_litter_N[days_gone+day]       = 0.0;
    all_out.mycorrhiza_litter_N[days_gone+day] = 0.0;
  }

  // leaf_mass_correction: leaf_mass becomes a genuinely accumulating pool
  // (growth.cpp's tree_state.needles in - the real, Makela-reduced and
  // sugar/nitrogen-gated needle growth demand, see that file's comment -
  // senescence-driven litter out) instead of fAPAR_and_LAI.cpp's
  // memoryless, annually-reset formula. Litter timing is a skew-normal
  // density over day-of-year (2026-11-10, replacing the earlier
  // fS^2-accumulator-plus-Makela-bell design) - see
  // development_notebooks/leaf_mass_correction.Rmd Sect. 5.1 for the
  // comparison behind this choice (validation RMSE 0.068 vs. 0.153 for
  // the old mechanism) and needle_litter_fS2_implementation_plan.md for
  // why that older mechanism itself existed (fixed a 124-day
  // dormant-season smear and mid-summer-wobble false triggers in an even
  // earlier, pre-fS^2 formula).
  if (boolsettings.leaf_mass_correction) {
    // day is 0-indexed here (day==0 is Jan 1); day_of_year is converted
    // to the 1-indexed convention (1-366) the R-side fit used.
    double day_of_year = static_cast<double>(day) + 1.0;

    // Detect the first day this year fS drops below 0.5 (searched only
    // after day_of_year 150, to skip the mid-summer wobble - fS sitting
    // at ~0.83-0.87 for 2-3 weeks some years - which a threshold close to
    // 1.0 would catch instead of the real autumn decline; only fS<0.5
    // reliably lands on it, confirmed in the Sect. 5.1 comparison).
    // Frozen at day==0 (growth.cpp) into fS_cross_day_last_year/
    // xi_litter for NEXT year's use - this year's own crossing can't be
    // known yet this early in a forward, day-by-day simulation.
    if (!tree_state.fS_crossed_this_year && day_of_year > 150 && photosynthesis.fS < 0.5) {
      tree_state.fS_crossed_this_year = true;
      tree_state.fS_cross_day_this_year = day_of_year;
    }

    // Skew-normal density: xi_litter (growth.cpp, from last year's
    // crossing day + parameters.xi_offset_litter), omega_litter/
    // alpha_litter fixed (fitted once against 2006-2016 Hyytiala
    // litter-trap data, see the Rmd Sect. 5.1). xi_litter stays -1 until
    // a prior year's crossing exists (bootstrap first year), giving 0
    // litter that year - same behaviour the sum_degree_fL_last_year<=1e-6
    // guard below already gave the old mechanism.
    double rate_unnormalized = 0.0;
    if (tree_state.xi_litter > 0) {
      double z = (day_of_year - tree_state.xi_litter) / parameters.omega_litter;
      double skew_density = (2.0 / parameters.omega_litter)
                             * (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * z * z)
                             * 0.5 * (1.0 + std::erf(parameters.alpha_litter * z / std::sqrt(2.0)));
      rate_unnormalized = skew_density;
    }
    tree_state.sum_degree_fL_running += rate_unnormalized;

    // One-year-lagged self-normalization - unchanged from the old
    // mechanism (function_structures.h). First simulated year: no
    // prior-year total yet, hold at 0.
    //
    // Magnitude (2026-11-10): the oldest needle cohort's real, tracked
    // mass (growth.cpp - the actual sugar/N-gated growth from n_age years
    // ago, not the N_repo allometric target), frozen at this year's
    // day==0 into needle_litter_magnitude_this_year. Replaces the old
    // N_repo_current/n_age * litter_boost_today correction entirely - the
    // cohort mass is already a real, bounded quantity (it can't shed more
    // than was actually grown that year), so the over/undershoot
    // correction that target-based approach needed no longer applies.
    double needle_litter_today = 0.0;
    if (tree_state.sum_degree_fL_last_year > 1e-6) {
      needle_litter_today = tree_state.needle_litter_magnitude_this_year
                             * rate_unnormalized / tree_state.sum_degree_fL_last_year;
    }

    all_out.culm_growth.leaf_mass[days_gone+day] = std::max(0.0, all_out.culm_growth.leaf_mass[index_ref] - needle_litter_today + tree_state.needles);

    all_out.culm_growth.needle_litter[days_gone+day] = all_out.culm_growth.needle_litter[index_ref] + needle_litter_today;
    all_out.needle_litter[days_gone+day] = needle_litter_today;

    double needle_litter_N_gross = parameters.theta_N_needles / parameters.carbon_share * needle_litter_today;
    needle_litter_N_retrans = needle_litter_N_gross * parameters.retrans_needles;
    all_out.needle_litter_N[days_gone+day] = needle_litter_N_gross - needle_litter_N_retrans;

    // Frozen once this year's needle growth has definitively ended
    // (growth.cpp's needle_growth_season_ended latch) - litter's own
    // correction (above) should react to a settled peak, not one that
    // could still be rising from growth/litter overlap near litter's
    // day>=212 gate (confirmed real in ~1/4 of simulated years).
    if (!tree_state.needle_growth_season_ended) {
      tree_state.max_leaf_mass_running = std::max(tree_state.max_leaf_mass_running, all_out.culm_growth.leaf_mass[days_gone+day]);
    }

    // Cone litter (masting cycle, growth.cpp's day==0 block): this
    // year's whole cone_litter_magnitude_this_year (0 except masting
    // years) spread evenly over a fixed day 121-273 window, not a fitted
    // continuous shape - RMSE-fitting a cone window was found to be
    // meaningless for a masting process (leaf_mass_correction.Rmd Sect.
    // 5.2: it just inflates duration chasing masting noise), so this
    // reuses that section's own fixed-window choice directly.
    double cone_litter_today = 0.0;
    if (day_of_year >= 121.0 && day_of_year <= 273.0) {
      cone_litter_today = tree_state.cone_litter_magnitude_this_year / 153.0;
    }
    all_out.cone_litter[days_gone+day] = cone_litter_today;
    all_out.cone_litter_N[days_gone+day] = 0.0;  // no N model for cones yet
    all_out.culm_growth.cone_litter[days_gone+day] = all_out.culm_growth.cone_litter[index_ref] + cone_litter_today;
  }

  // Credit for sugar_model.cpp's nitrogen_balance update, read on the
  // NEXT day's call (same one-day lag uptake.total_uptake already uses -
  // this function runs after sugar_model() each day, so today's litter
  // can't be credited into today's own balance). Always on, like
  // needle_litter itself - retranslocation isn't an alternate/
  // experimental mechanism to toggle, it's just how real litter N
  // accounting works.
  tree_state.N_retranslocated_today =
    needle_litter_N_retrans + root_litter_N_retrans + myco_litter_N_retrans;

  // Running within-year max of leaf_mass/roots, for root_leaf_correction's
  // yearly C_R and root_turnover's daily litter boost (see growth.cpp) -
  // tracked unconditionally since it's negligible cost and keeps it ready
  // if either setting is toggled on. max_leaf_mass_running_for_roots is a
  // dedicated second tracker, not a reuse of needles' own
  // max_leaf_mass_running - see that field's comment
  // (function_structures.h) for why reuse would read as already-reset 0.
  tree_state.max_leaf_mass_running_for_roots = std::max(tree_state.max_leaf_mass_running_for_roots, all_out.culm_growth.leaf_mass[days_gone+day]);
  tree_state.max_roots_running               = std::max(tree_state.max_roots_running,               all_out.culm_growth.roots[days_gone+day]);

  // Needles accumulation (potentially include drop logic later)

  if (day == 0) {
    all_out.culm_growth.needles[days_gone+day] = 0.0;
  } else {
    all_out.culm_growth.needles[days_gone+day] = all_out.culm_growth.needles[index_ref] + tree_state.needles;
  }

  // Bud accumulation - same lifetime-running-total pattern as needles
  // above, of the day's actual (resource-gated, line 58) bud growth.
  // culm_growth.bud was declared but never written anywhere until now.
  if (day == 0) {
    all_out.culm_growth.bud[days_gone+day] = 0.0;
  } else {
    all_out.culm_growth.bud[days_gone+day] = all_out.culm_growth.bud[index_ref] + tree_state.bud;
  }

  /*
   * Ring width and log
   */

  ring_width_generator(day,
                       days_gone,
                       tree_state,
                       all_out,
                       parameters);

}

