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
