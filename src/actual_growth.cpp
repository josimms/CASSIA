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

  double phloem_share   =   all_out.culm_growth.phloem[day + days_gone] / (all_out.culm_growth.phloem[day + days_gone] + all_out.culm_growth.xylem_sh[day + days_gone] + all_out.culm_growth.xylem_st[day + days_gone]);
  double xylem_st_share = all_out.culm_growth.xylem_st[day + days_gone] / (all_out.culm_growth.phloem[day + days_gone] + all_out.culm_growth.xylem_sh[day + days_gone] + all_out.culm_growth.xylem_st[day + days_gone]);
  double xylem_sh_share = all_out.culm_growth.xylem_sh[day + days_gone] / (all_out.culm_growth.phloem[day + days_gone] + all_out.culm_growth.xylem_sh[day + days_gone] + all_out.culm_growth.xylem_st[day + days_gone]);

  /*
   * Height
   */

  double storage_height = 0;
  if (boolsettings.sperling_model) {
    storage_height = (phloem_share * storage.phloem + xylem_st_share * storage.xylem_st + xylem_sh_share * storage.xylem_sh);
  } else {
    storage_height = storage.needles;
  }
  tree_state.height = tree_state.height * std::min(storage_height, nitrogen_capacity.height);

  /*
   * Wall
   */

  double storage_wall = 0;
  if (boolsettings.sperling_model) {
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
  if (boolsettings.sperling_model) {
    storage_roots = storage.roots;
  } else {
    storage_roots = storage.needles;
  }
  tree_state.roots = tree_state.roots * std::min(storage_roots, nitrogen_capacity.roots);

  /*
   * Mycorrhiza
   *
   * double storage_mycorrhiza;
   if (sperling_sugar_model) {
   storage_mycorrhiza = storage.roots;
   } else {
   storage_mycorrhiza = storage.needles;
   }
   actual_growth_out. = tree_state.roots * std::min(storage_mycorrhiza, nitrogen_capacity.roots);
   */

  /*
   * GD
   */
  double storage_GD;
  if (boolsettings.sperling_model) {
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

  all_out.culm_growth.roots[days_gone+day]      = all_out.culm_growth.roots[index_ref];
  all_out.culm_growth.mycorrhiza[days_gone+day] = all_out.culm_growth.mycorrhiza[index_ref];

  // Needles accumulation (potentially include drop logic later)

  if (day == 0) {
    all_out.culm_growth.needles[days_gone+day] = 0.0;
  } else {
    all_out.culm_growth.needles[days_gone+day] = all_out.culm_growth.needles[index_ref] + tree_state.needles;
  }


  std::cout << " tree_state.needles " << tree_state.needles;

  /*
   * Ring width and log
   */

  ring_width_generator(day,
                       days_gone,
                       tree_state,
                       all_out,
                       parameters);

}

