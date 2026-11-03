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

    double root_litter_today = all_out.culm_growth.roots[index_ref] * daily_turnover_fraction;
    double myco_litter_today = all_out.culm_growth.mycorrhiza[index_ref] * daily_turnover_fraction;

    all_out.culm_growth.roots[days_gone+day]      = std::max(0.0, all_out.culm_growth.roots[index_ref]      - root_litter_today + tree_state.roots);
    all_out.culm_growth.mycorrhiza[days_gone+day] = std::max(0.0, all_out.culm_growth.mycorrhiza[index_ref] - myco_litter_today + tree_state.ecto);

    all_out.culm_growth.root_litter[days_gone+day]       = all_out.culm_growth.root_litter[index_ref]       + root_litter_today;
    all_out.culm_growth.mycorrhiza_litter[days_gone+day] = all_out.culm_growth.mycorrhiza_litter[index_ref] + myco_litter_today;
  } else {
    all_out.culm_growth.roots[days_gone+day]      = all_out.culm_growth.roots[index_ref]      + tree_state.roots;
    all_out.culm_growth.mycorrhiza[days_gone+day] = all_out.culm_growth.mycorrhiza[index_ref] + tree_state.ecto;
  }

  // Needles accumulation (potentially include drop logic later)

  if (day == 0) {
    all_out.culm_growth.needles[days_gone+day] = 0.0;
  } else {
    all_out.culm_growth.needles[days_gone+day] = all_out.culm_growth.needles[index_ref] + tree_state.needles;
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

