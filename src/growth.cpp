#include "CASSIA.h"

/*
 * Growth Function
 */

void growth(int day,
            int days_gone,
            int year,

            growth_state& state,
            output_vector& out,

            double TAir,
            double TSoil_A,
            double TSoil_B,
            double Soil_Moisture,
            double PF,
            double GPP_ref,
            Settings boolsettings,

            const CASSIA_common& common,
            const CASSIA_parameters& parameters,
            const CASSIA_ratios& ratio,

            double CH,
            double B0,
            double GPP_mean,
            double GPP_previous_sum,

            double last_year_HH,
            int no_day)
{

  if (TAir > 0) {
    /*
     * Growth
     */
    state.g = 1 / (1 + exp(-common.a * (TAir - common.b)));
    if (state.g < 0) {
      state.g = 0;
    }
  }

  /*
   * Height
   */
  if (day == 0) {
    state.sH = parameters.sH0 + state.g; // First day that it could have a value is this!
  } else {
    state.sH = state.sH + state.g;
  }
  if ((state.sH > 0) & (state.sH < parameters.sHc)) {
    state.fH = (sin(2 * M_PI/parameters.sHc * (state.sH - parameters.sHc / 4)) + 1) / 2;
  } else {
    state.fH = 0;
  }
  double LH = parameters.LH0 * ratio.height_growth_coefficient;
  if (boolsettings.LH_estim) {
    LH = LH * GPP_previous_sum / parameters.GPP_mean;
  }
  if (boolsettings.tests) {
    if ( year== 2016) {
      LH = 37.70821;
    } else if (year == 2017) {
      LH = 35.35192;
    } else if (year == 2018) {
      LH = 36.93803;
    }
  }
  state.GH = state.g * state.fH * LH;
  state.height = 0.02405282 * 200.0 * state.GH / 1000.0 * ratio.form_factor;

  if (day == 0) {
    state.HH = parameters.HH0 + state.GH;
  } else {
    state.HH = state.HH + state.GH;
  }

  /*
   * Needles
   */
  double old_GN = 0.0;
  if (day == 0) {
    state.sN = parameters.sN0 + state.g;
  } else {
    state.sN = state.sN + state.g;
  }
  if ((state.sN > 0) & (state.sN < pow(parameters.sNc, 2))) {
    state.fN = (parameters.sNc * pow(state.sN, 0.5) - state.sN) / (pow(parameters.sNc, 2.0)/4.0);
  } else {
    state.fN = 0.0;
  }
  double LN = parameters.LN0;
  if (boolsettings.LN_estim) {
    LN = parameters.LN0 * GPP_previous_sum / parameters.GPP_mean;
  }
  if (boolsettings.tests) {
    if ( year== 2016) {
      LN = 1.971561;
    } else if (year == 2017) {
      LN = 1.848363;
    } else if (year == 2018) {
      LN = 1.931293;
    }
  }
  state.GN = state.g * state.fN * LN;

  double cumsum_GN = state.GN + old_GN; // TODO: should be initalised properly, but currently really really small values
  double HN = parameters.HN0 + cumsum_GN; // TODO: where does this go?

  state.needles = common.m_N * state.GN * last_year_HH / parameters.h_increment;

  state.max_N = std::max(state.needles, state.max_N);

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
  state.bud = state.g * fB * parameters.LB;

  /*
   * Diameter
   */

  if (day == 0) {
    state.sD = parameters.sD0_Trad;
  } else if (day <= parameters.diameter_start_day) {
    state.sD = state.sD;
  } else if (day >= parameters.diameter_start_day+1) {
    if (boolsettings.sD_estim_T_count) {
      // TODO: not coded! need to add!
      state.sD = state.sD + state.g;
    } else {
      state.sD = state.sD + state.g;
    }
  }
  if ((state.sD > 0) & (state.sD < pow(parameters.sDc, 2.0))) {
    state.fD = (parameters.sDc*pow(state.sD, 0.5) - state.sD)/(pow(parameters.sDc, 2.0)/4);
  } else {
    state.fD = 0;
  }

  double LD = parameters.LD0 * ratio.diameter_growth_coefficient;
  if (boolsettings.LD_estim) {
    if (day == 0) {
      state.S_GPP = 0.0;
      state.S_GPP_ref = 0.0;
    } else {
      state.S_GPP = state.S_GPP + state.dS_GPP;
      state.S_GPP_ref = state.S_GPP_ref + state.dS_GPP_ref;
    }
    state.dS_GPP = (PF - state.S_GPP) / parameters.tau_GPP;
    state.dS_GPP_ref = (GPP_ref - state.S_GPP_ref) / parameters.tau_GPP;

    // Daily LD depends on the GPP of five previous days:
    if (day > 78) { // TODO: hard coded! Should this be parameter diameter start day?
      LD = parameters.LD0 * ratio.diameter_growth_coefficient * state.S_GPP / state.S_GPP_ref;
    }
  }

  double old_GD = state.GD;
  state.GD = state.g * state.fD * LD;
  double tot_cells_pot = old_GD + state.GD;

  // "Uggla" divides cells to early and late wood
  if (state.sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    state.tau_E = parameters.tau_Ee;
    state.tau_W = parameters.tau_We;
  } else {
    state.tau_E = parameters.tau_El;
    state.tau_W = parameters.tau_Wl;
  }

  // Number of cells in enlargement, wall formation and mature phases on each day
  // TODO: why doesn't this section use the updated variables?
  double n_E_pot{state.n_E_pot}, n_W_pot{state.n_W_pot}, n_M_pot{state.n_M_pot};
  if (day < 1) {
    state.n_E_pot = 0.0;
    state.n_W_pot = 0.0;
    state.n_M_pot = 0.0;
  } else {
    state.n_E_pot = n_E_pot + state.GD - n_E_pot / state.tau_E;
    state.n_W_pot = n_W_pot + n_E_pot / state.tau_E - n_W_pot / state.tau_W;
    state.n_M_pot = n_M_pot + n_W_pot / state.tau_W;
  }

  // Carbon to enlargement of one earlywood/latewood cell per one day and cell (kg C cell-1 day-1)
  double CE_ew = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_ew / 2, 2) * 0.00266 * 342.296 / (common.gas_const * (TAir + common.abs_zero) * 1000) * 12 * common.M_C / 342.296 / state.tau_E;
  // TODO: should it be tau_W
  double CE_lw = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_lw / 2, 2) * 0.00266 * 342.296 / (common.gas_const * (TAir + common.abs_zero) * 1000) * 12 * common.M_C / 342.296 / state.tau_E;
  // The number of forming cell rows in the tree
  state.n_rows = ratio.form_factor * parameters.h0 / 0.00266 * M_PI * 0.175 / parameters.cell_d_ew;
  double carbon_enlargement_pot;
  if (state.sD < pow(parameters.sDc, 2) / common.Uggla) {
    carbon_enlargement_pot = CE_ew * state.n_E_pot;
  } else {
    carbon_enlargement_pot = CE_lw * state.n_E_pot;
  }

  double tau_Ee_floor = std::ceil(parameters.tau_Ee);
  state.use = state.n_rows * carbon_enlargement_pot;
  state.release = 0.0;

  // Carbon to wall formation
  // Carbon.daily.rate determined in parameters_common.R but NOTE!!!! not used at the moment, replaced by a parameter set to result in density app. 200 kg C m-3!
  double CW;
  if (state.sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    // CW = Carbon.daily.rate.ew;
    CW = 1.8e-11;
  } else {
    // CW = Carbon.daily.rate.lw;
    CW = 1.8e-11;
  }

  // The use of carbon to wall growth kg C per day
  state.diameter = state.n_rows * CW * state.n_W_pot;
  double cells_pot = state.n_W_pot + state.n_M_pot;

  double ew_cells_pot, lw_cells_pot;
  double ew_cells_pot_max;
  if (state.sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    ew_cells_pot = cells_pot;
    ew_cells_pot_max = std::max(ew_cells_pot, state.ew_cells_pot_max);
  } else {
    ew_cells_pot = 0.0;
    ew_cells_pot_max = std::max(ew_cells_pot, state.ew_cells_pot_max);

  }
  if (state.sD > pow(parameters.sDc, 2.0) / common.Uggla) {
    lw_cells_pot = cells_pot - ew_cells_pot - state.ew_cells_pot_max;
  } else {
    lw_cells_pot = 0.0;
  }
  double ew_width_pot = ew_cells_pot * parameters.cell_d_ew * 1000;
  double lw_width_pot = lw_cells_pot * parameters.cell_d_lw * 1000;

  double ew_mass_pot = ew_cells_pot * parameters.cell_volume_growth_per_day_ew * parameters.cell_wall_density_ew * 1000;
  double lw_mass_pot = lw_cells_pot * parameters.cell_volume_growth_per_day_lw * parameters.cell_wall_density_lw * 1000;

  if (state.sD < pow(parameters.sDc, 2.0) / common.Uggla) {
    state.pot_mm = ew_width_pot;
    state.pot_mm_max = std::max(state.pot_mm, state.pot_mm_max);
    // pot_mass_max = pot_mm_max * parameters.cell_volume_growth_per_day_ew * parameters.cell_wall_density_ew * 1000;
  } else {
    state.pot_mm_max = std::max(state.pot_mm, state.pot_mm_max);
    state.pot_mm = lw_width_pot + state.pot_mm_max;
    // pot_mass_max = lw_width_pot * parameters.cell_volume_growth_per_day_lw * parameters.cell_wall_density_lw * 1000; // TODO: fix this!
  }

  /*
   * Roots
   */
  double gR_fib, gR_pio, gR, LR;

  if (boolsettings.root_as_Ding) {
    double fib_coef = 0.25;                                  // 0.25, 2.3  # Determines the proportion of fibrous roots (1 leads to 37 % of fibrous roots, 0.25 to 13 % of fibrous roots and 2.3 to 63 % of fibrous roots)
    if ((day >= 149) & (day < 319)) {
      state.fR = 1.0/(1.0+exp(-0.038014*(day-148-56.06243)));
    } else {
      state.fR = 0.0;
    }

    gR_fib = -0.84 + 0.13 * TSoil_A - 0.44 + 2.11 * Soil_Moisture;              // growth of fibrous roots from Ding et al. 2019 (model 5)
    if (gR_fib < 0) {gR_fib = 0;}
    gR_pio = -0.84 + 0.13 * TSoil_B + 0.32 - 0.16 + 0.78 * Soil_Moisture;       // growth of pioneer roots from Ding et al. 2019 (model 5)
    if (gR_pio < 0) {gR_pio = 0;}
    gR = fib_coef * gR_fib + gR_pio;                                    // if fib_coef = 1 this leads in year 2018 to 37 % fibrous roots of all roots
    LR = 0.0049 * 1.0 / (0.37 * fib_coef + 0.63);                  // if fib_coef = 1 this leads to (roughly and on average) same total root growth as original. If fib_coef is changed, L is changed accordingly
  } else {
    LR = parameters.LR0 / parameters.root_lifetime;
    if (TSoil_B > common.TR0) {
      gR = (1 / (1 + exp(-common.a * ((TSoil_B - common.TR0) - common.b)))) * (1 - 1 / exp(Soil_Moisture * 10));
    } else {
      gR = 0;
    }

    if (day == 0) {
      state.sR = parameters.sR0 + state.g;
    } else {
      state.sR = state.sR + state.g;
    }
    if ((state.sR > 0) & (state.sR < parameters.sRc)) {
      state.fR = (sin(2 * M_PI / parameters.sRc * (state.sR - parameters.sRc / 4)) + 1) / 2;
    }
  }

  double GR = state.fR * LR * gR;
  if (GR > 0) {
    state.roots = GR;
  } else {
    state.roots = 0;
  }
  state.ecto = GR;

  /*
   * Save vector output
   */

  log_potential_growth(day, days_gone, state, out);
}

