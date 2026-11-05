#include "CASSIA.h"

// Makela (1999) eqn 28: C_j = [X_bal/X_current]^tau, capped at 1. Every
// other growth-rate factor in CASSIA (k(t), g_i(t), f_i(s_i(t)) in
// Schiestl-Aalto et al. 2015 eqn 4; the sugar/nitrogen capacity gates
// elsewhere in this codebase) only ever reduces growth below its
// calibrated maximum rate (L_i), never boosts it above - capping here
// keeps this correction damping-only, consistent with that architecture,
// rather than letting an undergrown organ exceed its calibrated ceiling.
double makela_control_factor(double X_bal, double X_current, double tau) {
  const double epsilon = 1e-6;
  if (X_current < epsilon || X_bal < epsilon) {
    return 1.0;
  }
  return std::min(1.0, std::pow(X_bal / X_current, tau));
}

// Mirror image of makela_control_factor, for LOSS rather than growth -
// but deliberately NOT capped the same damping-only way. Growth
// (makela_control_factor) stays damping-only: it never exceeds the
// calibrated baseline rate (1/n_age of N_repo/year), by design, same as
// every other growth-rate factor in this codebase (see that function's
// comment). Litter instead needs to move BOTH ways, by the reciprocal
// ratio: BOOSTS (> 1) when X_current (the pool) has overgrown past
// X_bal (N_repo) - shed faster, not just grow slower - and REDUCES
// (< 1) when the pool is undershooting N_repo, so an already-undergrown
// pool doesn't keep shedding at the full baseline rate while growth
// (capped, unable to compensate by growing faster) can never catch up.
// Without this, growth and litter both saturate at exactly the
// baseline rate once the pool sits below N_repo, locking the whole
// system into a stable cycle around N_repo/n_age that can never track a
// N_repo target that keeps rising as the tree grows (confirmed via
// debug trace, 1998-2022: peak tracks N_repo/n_age's trend (+0.025/yr),
// not N_repo's own (+0.076/yr) - a gap that widens from -0.12 to +1.18
// over 25 years). No floor/ceiling here at all - X_current/X_bal raised
// to tau is used directly, so the correction is genuinely symmetric.
double makela_loss_boost_factor(double X_bal, double X_current, double tau) {
  const double epsilon = 1e-6;
  if (X_current < epsilon || X_bal < epsilon) {
    return 1.0;
  }
  return std::pow(X_current / X_bal, tau);
}

// Same functional form as fN (growth.cpp's needle phenology curve) -
// zero at s=0, peaks at 1 when s=sLc^2/4, returns to zero past s=sLc^2.
// Self-terminating by construction, no external day-count cutoff needed.
double makela_bell_shape(double s, double sLc) {
  if (s <= 0 || s >= sLc * sLc) return 0.0;
  return (sLc * std::sqrt(s) - s) / (sLc * sLc / 4.0);
}

// Tested uncapped (2026-11-04) and reverted: without makela_control_factor's
// own cap, C_N_year/litter_boost_year together swing the pool between
// exactly 0 and ~1e15 within a few years. That cap on GROWTH stays
// (damping-only, as above). litter_boost_year's own cap was reinstated
// then too - but removing just that cap was retried on 2026-11-05 (see
// that function's comment) and is unstable for a DIFFERENT, more subtle
// reason than the blowup above: computed once per year (day == 0) from
// LAST year's peak, it's a one-year-LAGGED signal - a mild overshoot one
// year triggers a large litter cut the next, which undershoots, which
// triggers almost no litter the year after, which overshoots again -
// confirmed oscillating (pool swinging ~0 to ~7, N_repo ~6) rather than
// converging. The fix actually adopted: keep litter_boost_year's ratio
// uncapped/symmetric (growth.cpp's makela_loss_boost_factor), but stop
// computing it here from last year's frozen peak - actual_growth.cpp
// computes it fresh each day from THIS year's own max_leaf_mass_running
// (updated daily, reset at day == 0 below), so litter reacts to the
// season actually in progress, not a year-old snapshot. No lag, no
// oscillation.

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
            int no_day,
            double N_repo)
{

  if (TAir > 0) {
    /*
     * Growth
     */
    state.g = 1 / (1 + exp(-common.a * (TAir - common.b)));
    if (state.g < 0) {
      state.g = 0;
    }
  } else {
    state.g = 0;
  }

  // Makela (1999) structural-balance target: height's potential growth
  // corrects (damping only, see makela_control_factor above) toward
  // H_bal=alpha_H*D. Diameter is the uncorrected reference, not the
  // other way around: unlike height, diameter already has its own
  // within-season carbon-feedback (LD_estim/Eqn 19 of Schiestl-Aalto
  // et al. 2015, the daily current-year-GPP-vs-reference-year ratio
  // below) and the full cell-enlargement/wall-formation xylogenesis
  // cascade (ring_width.cpp) - a cruder external correction on top of
  // that more detailed, already-calibrated mechanism made diameter's
  // fit worse, not better (tested; see manuscript/paper_plan.Rmd
  // Working Notes). Needles/roots/bud were never corrected either (same
  // Working Notes: needle mass resets near-zero each season, blowing up
  // the control ratio; roots genuinely vary with nitrogen availability,
  // not a fixed ratio; bud has no analogous structural variable).
  // alpha_H = 89.45 (mean of two real SMEAR height/diameter series fits,
  // 88.8/90.1). tau_H comes from Settings (R-level, calibratable). Uses
  // yesterday's tracked state (today's culm_growth isn't written yet).
  double H_bal = 0.0, H_current = 0.0;
  if (boolsettings.height_diameter_correction) {
    const double alpha_H = 89.45;

    int idx_ref = days_gone + day - 1;
    if (idx_ref < 0) {
      idx_ref = 0;
    }
    double D_current = out.culm_growth.diameter[idx_ref];
    H_current = out.culm_growth.height[idx_ref];
    H_bal = alpha_H * D_current;
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
  if (boolsettings.height_diameter_correction) {
    LH = LH * makela_control_factor(H_bal, H_current, boolsettings.tau_H);
  }
  state.GH = state.g * state.fH * LH;
  state.height = B0 * CH * state.GH / 1000.0 * ratio.form_factor;

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
    state.needle_growth_season_ended = false;
  } else {
    state.sN = state.sN + state.g;
  }
  if ((state.sN > 0) & (state.sN < pow(parameters.sNc, 2))) {
    state.fN = (parameters.sNc * pow(state.sN, 0.5) - state.sN) / (pow(parameters.sNc, 2.0)/4.0);
  } else {
    state.fN = 0.0;
  }
  // One-way latch, true for the rest of the year once sN has crossed
  // sNc^2 - i.e. once this year's needle growth (fN) has definitively
  // finished, not just "hasn't started yet" (fN is also 0 before sN
  // turns positive in early spring - state.sN > 0's check above excludes
  // that case here). sN is monotonically non-decreasing within a year
  // (state.g >= 0 always), so this can only flip once per year. Used in
  // actual_growth.cpp to freeze max_leaf_mass_running once growth has
  // genuinely ended, rather than letting litter's correction react to a
  // reference that could still be rising from residual growth-season
  // overlap - confirmed real (2026-11-05 debug trace): fN is still
  // nonzero at day 212 (litter's own gate) in 6 of 25 simulated years.
  if (!state.needle_growth_season_ended && state.sN >= pow(parameters.sNc, 2)) {
    state.needle_growth_season_ended = true;
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

  // Makela (1999) structural-balance correction (2026-11-05 redesign),
  // applied DIRECTLY to state.needles - the real needle-growth demand
  // signal that sugar_model.cpp already uses in its carbon/nitrogen
  // accounting (growth_resp.needles, sugar_to_invest, nitrogen_target
  // all reference tree_state.needles). An earlier design computed a
  // SEPARATE, parallel growth term (state.leaf_mass_growth, scaled off
  // N_repo/fN_normalized instead of this real signal) - confirmed NOT
  // self-consistent with the carbon economy: sugar_model.cpp never
  // referenced leaf_mass_growth at all, so that pathway's needle growth
  // never competed for carbon with height/diameter/root growth, nor was
  // it ever actually charged against the sugar pool - exactly why
  // letting it boost (tried the same day) had no real cost anywhere in
  // the model. Reducing state.needles HERE instead, before sugar_model()
  // runs, makes the correction flow through the existing accounting
  // properly, and makes it damping-only by construction - needle growth
  // can only ever be reduced below what GN/fN/sugar availability would
  // otherwise produce, never inflated past it (makela_control_factor is
  // capped at 1, same as height's own correction and every other
  // growth-rate factor in this codebase).
  //
  // Accepted consequence (2026-11-05, on request): state.needles's own
  // sugar-gated annual total already runs well below N_repo/n_age in
  // this configuration (~0.3-0.8 vs ~1.4-2.0 kg/yr) - a damping-only
  // correction can't close that gap, only widen it further if the pool
  // ever overshoots. leaf_mass will likely persistently undershoot
  // N_repo under this design. That's an honest reflection of the real
  // carbon constraint already present in the model, not something to
  // paper over with an uncosted parallel growth term.
  //
  // C_N_year is computed ONCE per year, at day == 0, from LAST year's
  // settled peak - NOT recomputed daily from this year's own rising
  // peak (which is what this notebook tried first, same-day as litter's
  // own fix). Litter's correction has to be same-year/daily because it's
  // symmetric/uncapped (a lagged symmetric correction was confirmed
  // unstable, see actual_growth.cpp) - but growth's is capped at 1, so a
  // lagged value can't run away the same way, and a single yearly
  // constant avoids a real risk a daily one has: if the pool ever
  // overshot N_repo WHILE growth is still active (fN > 0, i.e. before
  // ~day 200), a same-year-reactive C_N_year would start damping growth
  // mid-season, continuously changing as the peak rises - distorting
  // fN's intended within-season shape (bounded, since capped, but real).
  // A yearly value multiplies that shape by one constant instead -
  // shape stays correct, only the year's overall magnitude changes. In
  // this particular run it makes no numerical difference either way
  // (C_N_year sits at exactly 1 all year regardless, confirmed via
  // isolation test - the pool never overshoots), but removes that latent
  // risk for conditions where it might.
  if (boolsettings.leaf_mass_correction && day == 0) {
    double leaf_mass_max_last_year = state.max_leaf_mass_running;
    state.C_N_year = makela_control_factor(N_repo, leaf_mass_max_last_year, boolsettings.tau_N);
    state.max_leaf_mass_running = 0.0;

    // Freeze/reset for needle litter's own one-year-lagged normalizer
    // (sum_degree_fL_running/last_year, function_structures.h).
    if (state.sum_degree_fL_running > 1e-6) {
      state.sum_degree_fL_last_year = state.sum_degree_fL_running;
    }
    state.sum_degree_fL_running = 0.0;
    state.sL_litter = 0.0;  // fresh accumulator for the new year
  }

  if (boolsettings.leaf_mass_correction) {
    state.needles *= state.C_N_year;
    state.N_repo_current = N_repo;
  }

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
  double CE_ew = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_ew / 2.0, 2.0) * 0.00266 * 342.296 / (common.gas_const * (TAir + common.abs_zero) * 1000.0) * 12.0 * common.M_C / 342.296 / state.tau_E;
  // TODO: should it be tau_W
  double CE_lw = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_lw / 2.0, 2.0) * 0.00266 * 342.296 / (common.gas_const * (TAir + common.abs_zero) * 1000.0) * 12.0 * common.M_C / 342.296 / state.tau_E;
  // The number of forming cell rows in the tree

  // TODO: are the cell.l parameters parameters in the original or are they these values
  // TODO: surely this D0 isn't correct
  state.n_rows = ratio.form_factor * parameters.h0 / 0.00266 * M_PI * 0.175 / parameters.cell_d_ew;
  double carbon_enlargement_pot;
  if (state.sD < pow(parameters.sDc, 2.0) / common.Uggla) {
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
  // Wall potential growth in original code
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

  // Makela (1999) structural-balance target for roots, referenced against
  // leaf area rather than diameter (see manuscript/paper_plan.Rmd):
  // leaf area already has genuine bounded within-year and across-year
  // turnover (Tian 2021 cohort aging/senescence in compute_fAPAR_used()),
  // unlike diameter, so it gives roots something real to correct toward.
  // Unlike height's correction, this is recomputed only once per year
  // (day == 0) and held fixed for the rest of that year via
  // state.C_R_year, rather than chasing leaf area's within-year seasonal
  // swings daily - closer to how LH_estim/LD_estim's own annual terms
  // behave than to height's own daily recorrection.
  //
  // Uses the PREVIOUS YEAR'S MAXIMUM leaf_area/roots, not its last-day
  // value: state.max_leaf_area_running/max_roots_running are updated
  // daily in actual_growth.cpp and only reset below, so at this point
  // (day == 0, before today's actual_growth() call) they still hold the
  // just-finished year's running max - its peak, not wherever the pool
  // happened to sit at the literal year-end snapshot (which can
  // undershoot the peak once senescence/litter are active).
  //
  // alpha_R is an UNCALIBRATED PLACEHOLDER - unlike alpha_H (fit against
  // two real independent SMEAR height/diameter series), no equivalent
  // root-biomass-vs-leaf-area observational series has been sourced yet.
  // Do not treat results with root_leaf_correction = TRUE as calibrated
  // until a real alpha_R replaces this.
  if (boolsettings.root_leaf_correction && day == 0) {
    const double alpha_R = 0.15; // TODO: placeholder, needs real calibration data

    double leaf_area_max_last_year = state.max_leaf_area_running;
    double roots_max_last_year = state.max_roots_running;
    double R_bal = alpha_R * leaf_area_max_last_year;
    state.C_R_year = makela_control_factor(R_bal, roots_max_last_year, boolsettings.tau_R);

    // Start tracking the new year's max fresh.
    state.max_leaf_area_running = 0.0;
    state.max_roots_running = 0.0;
  }

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
  if (boolsettings.root_leaf_correction) {
    GR = GR * state.C_R_year;
  }
  if (GR > 0) {
    state.roots = GR;
    state.ecto = GR;
  } else {
    state.roots = 0;
    state.ecto = 0;
  }

  /*
   * Save vector output
   */

  log_potential_growth(day, days_gone, state, out);
}

