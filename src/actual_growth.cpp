#include "CASSIA.h"

growth_out actual_growth(CASSIA_parameters parameters,
                         CASSIA_common common,
                         carbo_tracker storage,
                         growth_out potential_growth,
                         respiration_out resp,
                         bool sperling_sugar_model) {

  growth_out actual_growth_out;

  /*
   * Height
   */
  actual_growth_out.height = potential_growth.height * (0.082179938 * storage.phloem + 0.821799379 * storage.xylem_st + 0.096020683 * storage.xylem_sh);

  /*
   * Wall
   */
  actual_growth_out.wall = potential_growth.diameter * (0.082179938 * storage.phloem + 0.821799379 * storage.xylem_st + 0.096020683 * storage.xylem_sh);

  /*
   * Bud
   */
  actual_growth_out.bud = storage.needles * potential_growth.bud;

  /*
   * Needles
   */
  actual_growth_out.needles = storage.needles * potential_growth.needles;

  /*
   * Roots
   */
  actual_growth_out.roots = storage.roots * potential_growth.roots;

  /*
   * GD
   */
  // TODO: storage is needles as the sugar model isn't being used here, but it should be a general storage term when the sugar model is being used
  actual_growth_out.GD = storage.needles * potential_growth.GD;

  /*
   * Respiration Growth
   */
  actual_growth_out.respiration_growth = common.Rg_R * actual_growth_out.roots +
    common.Rg_N * actual_growth_out.needles +
    common.Rg_S * (actual_growth_out.wall + actual_growth_out.height);

  /*
   * Respiration Maintenance
   */
  if (sperling_sugar_model) {
    double storage_needles_rm, storage_roots_rm, storage_xylem_sh_rm, storage_xylem_st_rm, storage_phloem_rm;
    if (storage.needles < 0.1) {
      storage_needles_rm = 0;
    } else {
      storage_needles_rm = 1;
    }

    if (storage.roots < 0.1) {
      storage_roots_rm = 0;
    } else {
      storage_roots_rm = 1;
    }

    if (storage.xylem_sh < 0.1) {
      storage_xylem_sh_rm = 0;
    } else {
      storage_xylem_sh_rm = 1;
    }

    if (storage.xylem_st < 0.1) {
      storage_xylem_st_rm = 0;
    } else {
      storage_xylem_st_rm = 1;
    }

    if (storage.phloem < 0.1) {
      storage_phloem_rm = 0;
    } else {
      storage_phloem_rm = 1;
    }

    // Respiration
    actual_growth_out.respiration_maintenance = storage_needles_rm * resp.RmN + storage.roots * resp.RmR + (0.082179938 * storage.phloem + 0.821799379 * storage.xylem_st + 0.096020683 * storage.xylem_sh) * resp.RmS;
  } else {
    actual_growth_out.respiration_maintenance = storage.respiration * resp.Rm_a;
  }

  return actual_growth_out;
};

