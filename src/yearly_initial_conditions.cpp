#include "CASSIA.h"

yearly_in yearly_initial_conditions(double days_per_year) {

  yearly_in out;

  /*
   * Growth initial conditions
   */

  out.g_sH.assign(0.0, days_per_year);
  out.sH.assign(0.0, days_per_year);
  out.fH.assign(0.0, days_per_year);
  out.GH.assign(0.0, days_per_year);
  out.cumsum_GH.assign(0.0, days_per_year);
  out.HH.assign(0.0, days_per_year);
  out.height_pot_growth.assign(0.0, days_per_year);
  out.LH.assign(0.0, days_per_year);

  out.g_sN_T.assign(0.0, days_per_year);
  out.sN.assign(0.0, days_per_year);
  out.fN.assign(0.0, days_per_year);
  out.GN.assign(0.0, days_per_year);
  out.cumsum_GN.assign(0.0, days_per_year);
  out.needle_pot_growth.assign(0.0, days_per_year);
  out.LN.assign(0.0, days_per_year);
  out.HN.assign(0.0, days_per_year);

  out.g_sD_T.assign(0.0, days_per_year);
  out.sD.assign(0.0, days_per_year);
  out.fD.assign(0.0, days_per_year);
  out.GD.assign(0.0, days_per_year);
  out.diameter_pot_growth.assign(0.0, days_per_year);
  out.LD.assign(0.0, days_per_year);

  out.GPP_ref.assign(0.0, days_per_year);
  out.S_GPP.assign(0.0, days_per_year);
  out.S_GPP_ref.assign(0.0, days_per_year);
  out.dS_GPP.assign(0.0, days_per_year);
  out.dS_GPP_ref.assign(0.0, days_per_year);

  out.g_R.assign(0.0, days_per_year);
  out.sR.assign(0.0, days_per_year);
  out.fR.assign(0.0, days_per_year);
  out.GR.assign(0.0, days_per_year);
  out.roots_pot_growth.assign(0.0, days_per_year);
  out.LR.assign(0.0, days_per_year);

  out.tot_cells_pot.assign(0.0, days_per_year);
  out.tot_cells.assign(0.0, days_per_year);
  out.n_rows.assign(0.0, days_per_year);

  out.g.assign(0.0, days_per_year);
  out.sRc.assign(0.0, days_per_year);

  /*
   * XYLOGENESIS INITIAL CONDITIONS
   */

  out.tau_E.assign(0.0, days_per_year);
  out.tau_W.assign(0.0, days_per_year);

  out.n_E.assign(0.0, days_per_year);
  out.n_W.assign(0.0, days_per_year);
  out.n_M.assign(0.0, days_per_year);

  out.n_E_pot.assign(0.0, days_per_year);
  out.n_W_pot.assign(0.0, days_per_year);
  out.n_M_pot.assign(0.0, days_per_year);

  out.GD.assign(0.0, days_per_year);

  out.CE_ew.assign(0.0, days_per_year);
  out.CE_lw.assign(0.0, days_per_year);
  out.CW.assign(0.0, days_per_year);

  out.carbon_enlargement.assign(0.0, days_per_year);
  out.wall_growth.assign(0.0, days_per_year);
  out.wall_tot.assign(0.0, days_per_year);

  out.ew_cells.assign(0.0, days_per_year);
  out.lw_cells.assign(0.0, days_per_year);

 return(out);
}


growth_values_out growth_values_out_init() {
  growth_values_out out;

  out.sH = 0.0;
  out.fH = 0.0;
  out.HH = 0.0;
  out.sN = 0.0;
  out.fN = 0.0;
  out.sD = 0.0;
  out.fD = 0.0;
  out.sR = 0.0;
  out.fR = 0.0;
  out.n_rows = 0.0;
  out.GH = 0.0;
  out.GN = 0.0;
  out.GD = 0.0;
  out.max_N = 0.0;
  out.S_GPP = 0.0;
  out.dS_GPP = 0.0;
  out.S_GPP_ref = 0.0;
  out.dS_GPP_ref = 0.0;
  out.GPP_ref = 0.0;
  out.ew_cells_pot_max = 0.0;
  out.en_pot_growth = 0.0;
  out.pot_mm_max = 0.0;
  out.wall_pot_growth = 0.0;
  out.n_E_pot = 0.0;
  out.n_W_pot = 0.0;
  out.n_M_pot = 0.0;
  out.n_E = 0.0;
  out.tau_E = 0.0;
  out.tau_W = 0.0;

  return(out);
}


carbo_tracker carbo_tracker_init() {
  carbo_tracker out;

  out.needles = 0.0; // Initial value from inputs
  out.phloem = 0.0; // Initial value from inputs
  out.xylem_sh = 0.0; // Initial value from inputs
  out.xylem_st = 0.0; // Initial value from inputs
  out.roots = 0.0; // Initial value from inputs
  out.mycorrhiza = 0.0; // Initial value 0
  out.respiration = 0.0;
  out.initial_amount = 0.0; // Initial value from inputs
  out.B = 0.0;

  return(out);
}

carbo_balance carbo_balance_init() {
  carbo_balance out;

  out.previous_values.Ad = carbo_tracker_init();
  out.previous_values.As = carbo_tracker_init();
  out.previous_values.storage_term = carbo_tracker_init();
  out.previous_values.sB0 = 0.0;
  out.previous_values.tree_alive = false;

  return(out);
}


ring_width_out ring_width_out_init() {
  ring_width_out out;

  out.n_E_tot = 0.0;
  out.n_W_tot = 0.0;
  out.n_M_tot = 0.0;
  out.ew_cells_tot = 0.0;
  out.tot_mm = 0.0;
  out.max_ew_cells_tot = 0.0;
  out.max_ew_width_tot = 0.0;

  return(out);
}
