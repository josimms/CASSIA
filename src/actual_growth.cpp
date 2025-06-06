#include "CASSIA.h"

growth_out actual_growth(CASSIA_parameters parameters,
                         CASSIA_common common,
                         carbo_tracker storage,
                         growth_out potential_growth,
                         respiration_out resp,
                         bool sperling_sugar_model,
                         growth_out nitrogen_capacity) {

  growth_out actual_growth_out;

  double phloem_mass = 7.410537931;
  double xylem_sh_mass = 74.10537931;
  double xylem_st_mass = 8.65862069;

  double phloem_respiration_share = phloem_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);
  double xylem_sh_respiration_share = xylem_sh_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);
  double xylem_st_respiration_share = xylem_st_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);

  /*
   * Height
   */

  double storage_height = 0;
  if (sperling_sugar_model) {
    storage_height = (phloem_respiration_share * storage.phloem + xylem_st_respiration_share * storage.xylem_st + xylem_sh_respiration_share * storage.xylem_sh);
  } else {
    storage_height = storage.needles;
  }
  actual_growth_out.height = potential_growth.height * std::min(storage_height, nitrogen_capacity.height);

  /*
   * Wall
   */

  double storage_wall = 0;
  if (sperling_sugar_model) {
    storage_wall = (phloem_respiration_share * storage.phloem + xylem_st_respiration_share * storage.xylem_st + xylem_sh_respiration_share * storage.xylem_sh);
  } else {
    storage_wall = storage.needles;
  }
  actual_growth_out.wall = potential_growth.diameter * std::min(storage_wall, nitrogen_capacity.wall);

  /*
   * Bud
   */
  actual_growth_out.bud = potential_growth.bud * std::min(storage.needles, nitrogen_capacity.bud);

  /*
   * Needles
   */
  actual_growth_out.needles = potential_growth.needles * std::min(storage.needles, nitrogen_capacity.needles);

  /*
   * Roots
   */
  double storage_roots;
  if (sperling_sugar_model) {
    storage_roots = storage.roots;
  } else {
    storage_roots = storage.needles;
  }
  actual_growth_out.roots = potential_growth.roots * std::min(storage_roots, nitrogen_capacity.roots);

  /*
   * Mycorrhiza
   *
   * double storage_mycorrhiza;
   if (sperling_sugar_model) {
   storage_mycorrhiza = storage.roots;
   } else {
   storage_mycorrhiza = storage.needles;
   }
   actual_growth_out. = potential_growth.roots * std::min(storage_mycorrhiza, nitrogen_capacity.roots);
   */

  /*
   * GD
   */
  double storage_GD;
  if (sperling_sugar_model) {
    storage_GD = (phloem_respiration_share * storage.phloem + xylem_st_respiration_share * storage.xylem_st + xylem_sh_respiration_share * storage.xylem_sh);
  } else {
    storage_GD = storage.needles;
  }
  actual_growth_out.GD = potential_growth.GD * std::min(storage_GD, nitrogen_capacity.wall);

  return actual_growth_out;
};

