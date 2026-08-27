#include "CASSIA.h"

repola_out repola(double diameter, double height, CASSIA_parameters parameters) {

  /*
   * parameters.b0_repo = -6.303;
  parameters.b1_repo = 14.472;
  parameters.b2_repo = -3.976;
  parameters.D0 = 0.175;
  parameters.n_age = 3;
   */

  double dski_repo = 2 + 1.25*diameter*100;

  // NOTE: needle mas could be an input
  double needle_mass = exp(parameters.b0_repo+parameters.b1_repo*dski_repo/(dski_repo+6) + parameters.b2_repo*height/(height+1));
  double m_N_tot = needle_mass*parameters.carbon_share;		// kg C / tree
  double m_N = m_N_tot/(parameters.n_age*parameters.n_length);				// youngest needles (kg C / mm needle)

  double m_R_tot = m_N_tot*0.5;				// kg C / tree // TODO: is this the roots?
  // VARIO changes from original code *0.7 # kg C / tree, keksitty!

  repola_out out;
  out.needle_mass = needle_mass;
  out.m_N_tot = m_N_tot;
  out.m_R_tot = m_R_tot;
  out.m_N = m_N;

  return out;
}
