#include "CASSIA.h"

respiration_out respiration(int day,
                            CASSIA_parameters parameters,
                            CASSIA_ratios ratios,
                            repola_out repola,
                            double TAir,
                            double TSoil,
                            bool temp_rise,
                            bool Rm_acclimation,
                            bool mN_varies,

                            // parameters that I am not sure about
                            double B0) {


  double m_S_tot = ratios.form_factor * B0 * parameters.h0 * parameters.density_tree * parameters.carbon_share;		// woody carbon mass
  double m_N_tot = 2.233819;
  double m_R_tot = 1.11691;
  double Ra_share = -0.0007 * pow(TSoil, 2.0) + 0.0424 * TSoil + 0.0273;					// share of autotrophic soil respiration
  if (Ra_share < 0) {
    Ra_share = 0;
  }

  double Rm_accl = 1;
  if ((temp_rise) & (Rm_acclimation)) {
    Rm_accl = 0.85;
  }

  double RmS_a = parameters.Rm_S * (std::exp(log(parameters.Q10_S) / 10 * TAir) - 0.7) * m_S_tot * ratios.sapwood_share * Rm_accl; // Maintenance respiration of wood
  double RmR_a = parameters.Rm_R * (std::exp(log(parameters.Q10_R) / 10 * TSoil) - std::exp(-log(parameters.Q10_R) / 2)) * m_R_tot * Ra_share * Rm_accl;	// Maintenance respiration of roots
  double RmN_a = parameters.Rm_N * (std::exp(log(parameters.Q10_N) / 10 * TAir) - 0.7) * m_N_tot * Rm_accl;						// Maintenance respiration of needles

  double m_N_tot2;
  if (mN_varies) {
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

  respiration_out out;
  out.RmN = RmN;
  out.RmS = RmS;
  out.RmR = RmR;
  out.Rm_a = RmN + RmS + RmR;

  return out;
}

/*
 * Test function
 */

// [[Rcpp::export]]
Rcpp::List respiration_test_cpp(Rcpp::DataFrame pCASSIA_parameters,
                                Rcpp::DataFrame pCASSIA_common,
                                Rcpp::DataFrame pCASSIA_ratios,
                                Rcpp::DataFrame pCASSIA_sperling,
                                std::vector<double> extras_sperling,
                                int day,
                                double TAir,
                                double TSoil,
                                bool temp_rise,
                                bool Rm_acclimation,
                                bool mN_varies,
                                double B0) {

  /*
   * Building the data structures
   */

  CASSIA_common common = make_common(pCASSIA_common);
  CASSIA_parameters parameters = make_CASSIA_parameters(pCASSIA_parameters, pCASSIA_sperling);
  CASSIA_ratios ratios = make_ratios(pCASSIA_ratios);
  repola_out repola_values = repola(parameters);

  respiration_out resp = respiration(day,
                                     parameters, ratios, repola_values,
                                     TAir, TSoil,
                                     temp_rise, Rm_acclimation, mN_varies,
                                     B0);

  return Rcpp::List::create(Rcpp::_["RmN"] = resp.RmN,
                            Rcpp::_["RmS"] = resp.RmS,
                            Rcpp::_["RmR"] = resp.RmR,
                            Rcpp::_["Rm_a"] = resp.Rm_a);
}


