#include "CASSIA.h"

ring_width_out ring_width_generator(int day,
                                    ring_width_out previous_value,
                                    growth_values_out growth_previous,
                                    CASSIA_parameters parameters,
                                    double GD_tot) {

  double n_E_tot, n_W_tot, n_M_tot;
  if (day < 1) {
    n_E_tot = 0;
    n_W_tot = 0;
    n_M_tot = 0;
  } else {
    n_E_tot = previous_value.n_E_tot + GD_tot - previous_value.n_E_tot / growth_previous.tau_E;
    n_W_tot = previous_value.n_W_tot + previous_value.n_E_tot / growth_previous.tau_E - previous_value.n_W_tot / growth_previous.tau_W;
    n_M_tot = previous_value.n_M_tot + previous_value.n_W_tot / growth_previous.tau_W;
  }

  double cells_tot = n_W_tot + n_M_tot;

  double ew_cells_tot, lw_cells_tot, max_ew_cells_tot, tot_mm;
  if (growth_previous.sD < pow(parameters.sDc, 2.0) / parameters.Uggla) {
    ew_cells_tot = cells_tot;
    lw_cells_tot = 0.0;
    max_ew_cells_tot = std::max(ew_cells_tot, previous_value.max_ew_cells_tot);
  } else {
    ew_cells_tot = 0.0;
    max_ew_cells_tot = std::max(ew_cells_tot, previous_value.max_ew_cells_tot);
    lw_cells_tot = cells_tot - ew_cells_tot - max_ew_cells_tot;
  }

  double ew_width_tot = ew_cells_tot * parameters.cell_d_ew * 1000;
  double lw_width_tot = lw_cells_tot * parameters.cell_d_lw * 1000;
  double max_ew_width_tot = std::max(ew_width_tot, previous_value.max_ew_width_tot);

  if (growth_previous.sD < pow(parameters.sDc, 2.0) / parameters.Uggla) {
    tot_mm = ew_width_tot;
  } else {
    tot_mm = lw_width_tot + max_ew_width_tot;
  }

  // Output
  ring_width_out out;
  out.n_E_tot = n_E_tot;
  out.n_M_tot = n_M_tot;
  out.n_W_tot = n_W_tot;
  out.tot_mm = tot_mm;
  out.max_ew_cells_tot = max_ew_cells_tot;
  out.max_ew_width_tot = max_ew_width_tot;

  return out;
}
