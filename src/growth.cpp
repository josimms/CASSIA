#include "CASSIA.h"

/*
 * Make Growth function
 */

growth_values_out make_growth(std::vector<double> values_in) {
  growth_values_out out;
  out.sH = values_in[0];
  out.fH = values_in[1];
  out.sN = values_in[2];
  out.fN = values_in[3];
  out.sD = values_in[4];
  out.fD = values_in[5];
  out.sR = values_in[6];
  out.fR = values_in[7];
  out.n_rows = values_in[8];
  out.GH = values_in[9];
  out.GN = values_in[10];
  out.S_GPP = values_in[11];
  out.dS_GPP = values_in[12];
  out.S_GPP_ref = values_in[13];
  out.dS_GPP_ref = values_in[14];
  out.GPP_ref = values_in[15];
  out.ew_cells_pot_max = values_in[16];
  out.pot_mm_max = values_in[17];
  out.wall_pot_growth = values_in[18];
  out.n_E_pot = values_in[19];
  out.n_W_pot = values_in[20];
  out.n_M_pot = values_in[21];
  out.n_E = values_in[22];

  return(out);
}


/*
 * Growth Function
 */

growth_out growth(int day,
                  int year,
                  double TAir,
                  double TSoil_A,
                  double TSoil_B,
                  double Soil_Moisture,
                  double PF,
                  double GPP_ref,
                  bool root_as_Ding,
                  bool xylogenesis_option,
                  bool environmental_effect_xylogenesis,
                  bool sD_estim_T_count,
                  CASSIA_common common,
                  CASSIA_parameters parameters,
                  CASSIA_ratios ratio,
                  double CH,
                  double B0,
                  double en_pot_growth_old,
                  double GPP_mean,
                  double GPP_previous_sum,

                  bool LH_estim,
                  bool LN_estim,
                  bool LD_estim,

                  growth_values_out growth_previous,
                  double last_year_HH,
                  int no_day)
{

  double g = 0;
  double g_sH = g;
  double g_sN = g;
  double g_sD_T;
  double sH, sN, sD, sR;
  double fH, fN, fD, fR;
  double GH, GN, GD, GR;
  double height_pot_growth, needle_pot_growth;
  double HH = 0;
  sH = sN = sD = sR = 0;
  fH = fN = fD = fR = 0;
  GH = GN = GD = GR = 0;
  height_pot_growth = needle_pot_growth = 0;

  if (TAir > 0) {
    /*
     * Growth
     */
    g = 1 / (1 + exp(-common.a * (TAir - common.b)));
    if (g < 0) {
      g = 0;
    }
  }

  /*
   * Height
   */
  g_sH = g;
  if (day == 0) {
    sH = parameters.sH0 + g_sH; // First day that it could have a value is this!
  } else {
    sH = growth_previous.sH + g_sH;
  }
  if (sH > 0 & sH < parameters.sHc) {
    fH = (sin(2 * M_PI/parameters.sHc * (sH - parameters.sHc / 4)) + 1) / 2;
  } else {
    fH = 0;
  }
  double LH = parameters.LH0 * ratio.height_growth_coefficient;
  if (LH_estim) {
    LH = LH * GPP_previous_sum / GPP_mean;
  }
  GH = g_sH * fH * LH;
  height_pot_growth = 0.02405282 * 200.0 * GH / 1000.0 * ratio.form_factor;

  if (day == 0) {
    HH = parameters.HH0 + GH;
  } else {
    HH = growth_previous.HH + GH;
  }

  /*
   * Needles
   */
  g_sN = g;
  if (day == 0) {
    sN = parameters.sN0 + g_sN;
  } else {
    sN = growth_previous.sN + g_sN;
  }
  if (sN > 0 & sN < pow(parameters.sNc, 2)) {
    fN = (parameters.sNc * pow(sN, 0.5) - sN) / (pow(parameters.sNc, 2.0)/4.0);
  } else {
    fN = 0.0;
  }
  double LN = parameters.LN0;
  if (LN_estim) {
    LN = parameters.LN0 * GPP_previous_sum / GPP_mean;
  }
  GN = g * fN * LN;

  double cumsum_GN = growth_previous.GN + GN;
  double HN = parameters.HN0 + cumsum_GN;

  needle_pot_growth = common.m_N * GN * last_year_HH / parameters.h_increment;

  double max_N = std::max(needle_pot_growth, growth_previous.max_N);

  /*
   * Bud growth
   */

  double fB = 0;
  double sB = day + 1 - parameters.sB0;
  if (day >= parameters.sB0) { // This value should come from the sugar model
    if (sB < parameters.sBc) {
      fB = (sin(2 * M_PI / parameters.sBc * (sB - parameters.sBc / 4)) + 1) / 2;
    }
  }
  double bud_pot_growth = g * fB * parameters.LB;

  /*
   * Diameter
   */

  g_sD_T = g;
  if (day == 0) {
    sD = parameters.sD0_Trad;
  } else if (day <= 78) {
    sD = growth_previous.sD;
  } else if (day >= 79) {
    if (sD_estim_T_count) {
      // TODO: not coded! need to add!
      sD = growth_previous.sD + g_sD_T;
    } else {
      sD = growth_previous.sD + g_sD_T;
    }
  }
  if (sD > 0 & sD < pow(parameters.sDc, 2.0)) {
    fD = (parameters.sDc*pow(sD, 0.5) - sD)/(pow(parameters.sDc, 2.0)/4);
  } else {
    fD = 0;
  }

  double LD = parameters.LD0 * ratio.diameter_growth_coefficient;
  double S_GPP, dS_GPP, S_GPP_ref, dS_GPP_ref;
  if (LD_estim) {
    LD = parameters.LD0 * ratio.diameter_growth_coefficient;
    if (day == 0) {
      S_GPP = 0.0;
      S_GPP_ref = 0.0;
    } else {
      S_GPP = growth_previous.S_GPP + growth_previous.dS_GPP;
      S_GPP_ref = growth_previous.S_GPP_ref + growth_previous.dS_GPP_ref;
    }
    dS_GPP = (PF - S_GPP) / parameters.tau_GPP;
    dS_GPP_ref = (GPP_ref - S_GPP_ref) / parameters.tau_GPP;

    // Daily LD depends on the GPP of five previous days:
    if (day > 78) {
      LD = parameters.LD0 * ratio.diameter_growth_coefficient * S_GPP / S_GPP_ref;
    }
  }

  GD = g_sD_T * fD * LD;
  double tot_cells_pot = GD + growth_previous.GD;

  double tau_E, tau_W;
  // "Uggla" divides cells to early and late wood
  if (sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    tau_E = parameters.tau_Ee;
    tau_W = parameters.tau_We;
  } else {
    tau_E = parameters.tau_El;
    tau_W = parameters.tau_Wl;
  }

  // Number of cells in enlargement, wall formation and mature phases on each day
  double n_E_pot, n_W_pot, n_M_pot;
  if (day < 1) {
    n_E_pot = 0.0;
    n_W_pot = 0.0;
    n_M_pot = 0.0;
  } else {
    n_E_pot = growth_previous.n_E_pot + GD - growth_previous.n_E_pot / tau_E;
    n_W_pot = growth_previous.n_W_pot + growth_previous.n_E_pot / tau_E - growth_previous.n_W_pot / growth_previous.tau_W;
    n_M_pot = growth_previous.n_M_pot + growth_previous.n_W_pot / growth_previous.tau_W;
  }

  // Carbon to enlargement of one earlywood/latewood cell per one day and cell (kg C cell-1 day-1)
  double CE_ew = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_ew / 2, 2) * 0.00266 * 342.296 / (common.gas_const * (TAir + common.abs_zero) * 1000) * 12 * common.M_C / 342.296 / tau_E;
  double CE_lw = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_lw / 2, 2) * 0.00266 * 342.296 / (common.gas_const * (TAir + common.abs_zero) * 1000) * 12 * common.M_C / 342.296 / tau_E;
  // The number of forming cell rows in the tree
  double n_rows = ratio.form_factor * parameters.h0 / 0.00266 * M_PI * 0.175 / parameters.cell_d_ew;
  double carbon_enlargement_pot;
  if (sD < pow(parameters.sDc, 2) / common.Uggla) {
    carbon_enlargement_pot = CE_ew * n_E_pot;
  } else {
    carbon_enlargement_pot = CE_lw * n_E_pot;
  }

  double tau_Ee_floor = std::ceil(parameters.tau_Ee);
  double en_pot_growth = n_rows * carbon_enlargement_pot;
  double en_pot_release = 0;
  if (day <= tau_Ee_floor) {
    en_pot_release = 0; // The carbon used in enlargement is released after some days.
  } else {
    en_pot_release = en_pot_growth_old;
  }

  // Carbon to wall formation
  // Carbon.daily.rate determined in parameters_common.R but NOTE!!!! not used at the moment, replaced by a parameter set to result in density app. 200 kg C m-3!
  double CW;
  if (sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    // CW = Carbon.daily.rate.ew;
    CW = 1.8e-11;
  } else {
    // CW = Carbon.daily.rate.lw;
    CW = 1.8e-11;
  }

  // The use of carbon to wall growth kg C per day
  double wall_pot_growth = n_rows * CW * n_W_pot;
  double cells_pot = n_W_pot + n_M_pot;

  double ew_cells_pot, lw_cells_pot;
  double ew_cells_pot_max;
  if (sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    ew_cells_pot = cells_pot;
    ew_cells_pot_max = std::max(ew_cells_pot, growth_previous.ew_cells_pot_max);
  } else {
    ew_cells_pot = 0.0;
    ew_cells_pot_max = std::max(ew_cells_pot, growth_previous.ew_cells_pot_max);

  }
  if (sD > pow(parameters.sDc, 2.0) / common.Uggla) {
    lw_cells_pot = cells_pot - ew_cells_pot - ew_cells_pot_max;
  } else {
    lw_cells_pot = 0.0;
  }
  double ew_width_pot = ew_cells_pot * parameters.cell_d_ew * 1000;
  double lw_width_pot = lw_cells_pot * parameters.cell_d_lw * 1000;
  double pot_mm, pot_mm_max;
  if (sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    pot_mm = ew_width_pot;
    pot_mm_max = std::max(pot_mm, growth_previous.pot_mm_max);
  } else {
    pot_mm_max = std::max(pot_mm, growth_previous.pot_mm_max);
    pot_mm = lw_width_pot + pot_mm_max;
  }

  /*
   * Roots
   */
  double gR_fib, gR_pio, gR, LR;

  if (root_as_Ding) {
    double fib_coef = 0.25;                                  // 0.25, 2.3  # Determines the proportion of fibrous roots (1 leads to 37 % of fibrous roots, 0.25 to 13 % of fibrous roots and 2.3 to 63 % of fibrous roots)
    if (day >= 149 & day < 319) {
      fR = 1/(1+exp(-0.038014*(day-148-56.06243)));
    } else {
      fR = 0;
    }

    gR_fib = -0.84 + 0.13 * TSoil_A -0.44 + 2.11 * Soil_Moisture;              // growth of fibrous roots from Ding et al. 2019 (model 5)
    if (gR_fib < 0) {gR_fib = 0;}
    gR_pio = -0.84 + 0.13 * TSoil_B + 0.32 -0.16 + 0.78 * Soil_Moisture;       // growth of pioneer roots from Ding et al. 2019 (model 5)
    if (gR_pio < 0) {gR_pio = 0;}
    gR = fib_coef * gR_fib + gR_pio;                                    // if fib_coef = 1 this leads in year 2018 to 37 % fibrous roots of all roots
    LR = 0.0049 * 1 / (0.37 * fib_coef + 0.63);                  // if fib_coef = 1 this leads to (roughly and on average) same total root growth as original. If fib_coef is changed, L is changed accordingly
  } else {
    LR = parameters.LR0 / parameters.root_lifetime;
    if (TSoil_B > common.TR0) {
      gR = (1 / (1 + exp(-common.a * ((TSoil_B - common.TR0) - common.b)))) * (1 - 1 / exp(Soil_Moisture * 10));
    } else {
      gR = 0;
    }

    double g_sR = g;
    if (day == 0) {
      sR = parameters.sR0 + g_sR;
    } else {
      sR = growth_previous.sR + g_sR;
    }
    if (sR > 0 & sR < parameters.sRc) {
      fR = (sin(2 * M_PI / parameters.sRc * (sR - parameters.sRc / 4)) + 1) / 2;
    }
  }

  GR = fR * LR * gR;
  double root_pot_growth;
  if (GR > 0) {
    root_pot_growth = GR;
  } else {
    root_pot_growth = 0;
  }

  /*
   * OUT!
   *
   * Creating the potential growth output
   */

  growth_values_out values_out;
  values_out.sH = sH;
  values_out.fH = fH;
  values_out.HH = HH;
  values_out.sN = sN;
  values_out.fN = fN;
  values_out.sD = sD;
  values_out.fD = fD;
  values_out.fR = fR;
  values_out.sR = sR;
  values_out.n_rows = n_rows;
  values_out.GH = GH;
  values_out.GN = GN;
  values_out.GD = GD;
  values_out.max_N = max_N;
  values_out.S_GPP = S_GPP;
  values_out.dS_GPP = dS_GPP;
  values_out.S_GPP_ref = S_GPP_ref;
  values_out.dS_GPP_ref = dS_GPP_ref;
  values_out.ew_cells_pot_max = ew_cells_pot_max;
  values_out.en_pot_growth = en_pot_growth;
  values_out.pot_mm_max = pot_mm_max;
  values_out.wall_pot_growth = wall_pot_growth;
  values_out.n_E_pot = n_E_pot;
  values_out.n_W_pot = n_W_pot;
  values_out.n_M_pot = n_M_pot;
  values_out.tau_E = tau_E;
  values_out.tau_W = tau_W;

  growth_out out;
  out.height = height_pot_growth;
  out.needles = needle_pot_growth;
  out.roots = root_pot_growth;
  out.diameter = wall_pot_growth;
  out.bud = bud_pot_growth;
  out.release = en_pot_release;
  out.use = en_pot_growth;
  out.previous_values = values_out;
  out.g = g;

  return out;
}


/*
 * Growth wrapper!
 */

// [[Rcpp::export]]
Rcpp::List growth_wrapper(int day,
                          int year,
                          double TAir,
                          double TSoil_A,
                          double TSoil_B,
                          double Soil_Moisture,
                          double PF,
                          double GPP_ref,
                          bool root_as_Ding,
                          bool xylogenesis_option,
                          bool environmental_effect_xylogenesis,
                          bool sD_estim_T_count,
                          Rcpp::DataFrame pCASSIA_common,
                          Rcpp::DataFrame pCASSIA_parameters,
                          Rcpp::DataFrame pCASSIA_ratios,
                          Rcpp::DataFrame pCASSIA_sperling,
                          std::vector<double> extras_sperling,

                          double CH,
                          double B0,
                          double en_pot_growth_old,
                          double GPP_mean,
                          double GPP_previous_sum,

                          bool LH_estim,
                          bool LN_estim,
                          bool LD_estim,

                          std::vector<double> growth_in,
                          double last_year_HH,
                          int no_day) {

  growth_values_out growth_previous = make_growth(growth_in);
  CASSIA_common common = make_common(pCASSIA_common);
  CASSIA_parameters parameters = make_CASSIA_parameters(pCASSIA_parameters, pCASSIA_sperling);
  CASSIA_ratios ratios = make_ratios(pCASSIA_ratios);

  growth_out out = growth(day,
                          year,
                          TAir,
                          TSoil_A,
                          TSoil_B,
                          Soil_Moisture,
                          PF,
                          GPP_ref,
                          root_as_Ding,
                          xylogenesis_option,
                          environmental_effect_xylogenesis,
                          sD_estim_T_count,
                          common,
                          parameters,
                          ratios,
                          CH,
                          B0,
                          en_pot_growth_old,
                          GPP_mean,
                          GPP_previous_sum,

                          LH_estim,
                          LN_estim,
                          LD_estim,

                          growth_previous,
                          last_year_HH,
                          no_day);

  return Rcpp::List::create(Rcpp::_["height_growth"] = out.height,
                            Rcpp::_["needles_growth"] = out.needles,
                            Rcpp::_["roots_growth"] = out.roots,
                            Rcpp::_["diameter"] = out.diameter,
                            Rcpp::_["bud"] = out.bud);

}




