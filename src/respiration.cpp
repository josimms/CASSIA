#include "CASSIA.h"

void respiration(int day,
                 growth_state& tree_state,
                 output_vector& out,
                 const CASSIA_parameters& parameters,
                 const CASSIA_ratios& ratios,
                 double TAir,
                 double TSoil,
                 const Settings& settings,
                 double B0) {


  double m_S_tot = ratios.form_factor * B0 * parameters.h0 * parameters.density_tree * parameters.carbon_share;		// woody carbon mass
  double m_N_tot = 2.233819;
  double m_R_tot = 1.11691;
  double Ra_share = -0.0007 * pow(TSoil, 2.0) + 0.0424 * TSoil + 0.0273;					// share of autotrophic soil respiration
  if (Ra_share < 0) {
    Ra_share = 0;
  }

  double Rm_accl = 1;
  if ((settings.temp_rise) & (settings.Rm_acclimation)) {
    Rm_accl = 0.85;
  }

  double RmS_a = parameters.Rm_S * (std::exp(log(parameters.Q10_S) / 10 * TAir) - 0.7) * m_S_tot * ratios.sapwood_share * Rm_accl; // Maintenance respiration of wood
  double RmR_a = parameters.Rm_R * (std::exp(log(parameters.Q10_R) / 10 * TSoil) - std::exp(-log(parameters.Q10_R) / 2)) * m_R_tot * Ra_share * Rm_accl;	// Maintenance respiration of roots
  double RmN_a = parameters.Rm_N * (std::exp(log(parameters.Q10_N) / 10 * TAir) - 0.7) * m_N_tot * Rm_accl;						// Maintenance respiration of needles

  double m_N_tot2;
  if (settings.mN_varies) {
    if ((day > 149) & (day < 284)) {
      m_N_tot2 = m_N_tot;
    } else {
      m_N_tot2 = m_N_tot * 2/3;
    }
    RmN_a = parameters.Rm_N * (std::exp(log(parameters.Q10_N) / 10 * TAir) - 0.7) * m_N_tot2 * Rm_accl;
  }

  double RmN = 0.0;
  double RmS = 0.0;
  double RmR = 0.0;
  if (RmN_a > 0.0) {
    RmN = RmN_a;
  }
  if (RmS_a > 0.0) {
    RmS = RmS_a;
  }
  if (RmR_a > 0.0) {
    RmR = RmR_a;
  }

  tree_state.RmN = RmN;
  tree_state.RmS = RmS;
  tree_state.RmR = RmR;
  tree_state.Rm_a = RmN + RmS + RmR;

  out.RmN.push_back(RmN);
  out.RmS.push_back(RmS);
  out.RmR.push_back(RmR);
  out.Rm_a.push_back(RmN + RmS + RmR);

}

