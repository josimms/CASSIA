#include "CASSIA.h"
#include <cmath>
#include <algorithm>

void ring_width_generator(int day,
                          int days_gone,
                          growth_state& state,
                          output_vector& all_out,
                          const CASSIA_parameters& parameters) {

  // double n_rows = ratios.form_factor * parameters.h0 / parameters.cell_l_ew * M_PI * parameters.D0 / parameters.cell_d_ew; // TODO: not used anywhere

  // Save previous values to avoid in-place dependency issues
  double prev_n_E = state.n_E_tot;
  double prev_n_W = state.n_W_tot;
  double prev_n_M = state.n_M_tot;
  double prev_max_ew_cells = state.max_ew_cells_tot;
  double prev_max_ew_width = state.max_ew_width_tot;

  if (day < 1) {
    state.n_E_tot = 0.0;
    state.n_W_tot = 0.0;
    state.n_M_tot = 0.0;
  } else {
    state.n_E_tot = prev_n_E + state.GD - prev_n_E / state.tau_E;
    state.n_W_tot = prev_n_W + prev_n_E / state.tau_E - prev_n_W / state.tau_W;
    state.n_M_tot = prev_n_M + prev_n_W / state.tau_W;
  }

  double cells_tot = state.n_W_tot + state.n_M_tot;

  if (state.sD < std::pow(parameters.sDc, 2.0) / parameters.Uggla) {
    state.ew_cells_tot = cells_tot;
    state.lw_cells_tot = 0.0;
    state.max_ew_cells_tot = std::max(state.ew_cells_tot, prev_max_ew_cells);
  } else {
    state.ew_cells_tot = 0.0;
    state.max_ew_cells_tot = std::max(state.ew_cells_tot, prev_max_ew_cells);
    state.lw_cells_tot = cells_tot - state.ew_cells_tot - state.max_ew_cells_tot;
  }

  double ew_width_tot = state.ew_cells_tot * parameters.cell_d_ew * 1000.0;
  double lw_width_tot = state.lw_cells_tot * parameters.cell_d_lw * 1000.0;
  state.max_ew_width_tot = std::max(ew_width_tot, prev_max_ew_width);

  if (state.sD < std::pow(parameters.sDc, 2.0) / parameters.Uggla) {
    state.tot_mm = ew_width_tot;
  } else {
    state.tot_mm = lw_width_tot + state.max_ew_width_tot;
  }

  /*
   * LOGGING
   */

  int index_ref = (day + days_gone) / 365;

  // NOTE: Already cumulative
  all_out.ring_width[day + days_gone] = state.tot_mm;

  all_out.culm_growth.diameter[day + days_gone] = all_out.culm_growth.diameter[index_ref] + 2.0*state.tot_mm/1000.0;

  all_out.culm_growth.diameter_potential[day + days_gone] = all_out.culm_growth.diameter_potential[index_ref] + 2.0*state.pot_mm/1000.0;

  /*
   * Commented out the ring width based calculations of the xylem and phloem growth!
   */

  // double xylem_width{0.0}, xylem_mass{0.0};
  // double days_in_16_years = 366*4 + 365*12;
  // if (day + days_gone > days_in_16_years) {
  // xylem_width = all_out.culm_growth.diameter[day + days_gone] - all_out.culm_growth.diameter[index_ref - days_in_16_years];
  // xylem_mass = M_PI * pow(xylem_width/100.0, 2.0) * all_out.culm_growth.height[index_ref] * parameters.cell_wall_density_ew;
  //   } else {
  // xylem_mass = 0.8 * M_PI * pow(all_out.ring_width[day + days_gone]/1000.0, 2.0) * all_out.culm_growth.height[index_ref] * parameters.cell_wall_density_ew;
  // }

  // double phloem_mass = M_PI * pow(1.5/1000.0, 2.0) * all_out.culm_growth.height[index_ref] * parameters.cell_wall_density_ew;

  // all_out.culm_growth.xylem_sh[day + days_gone] = all_out.culm_growth.xylem_sh[index_ref] + xylem_mass;
  // all_out.culm_growth.phloem[day + days_gone] = all_out.culm_growth.xylem_sh[index_ref] + phloem_mass;

}
