#include "CASSIA.h"

int leap_year(int year) {
  if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
    return 366;
  } else {
    return 365;
  }
};

#include <iostream>

bool is_leap_year(int year) {
  return (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
}

int days_in_month(int year, int month) {
  switch(month) {
  case 1: case 3: case 5: case 7: case 8: case 10: case 12:
    return 31;
  case 4: case 6: case 9: case 11:
    return 30;
  case 2:
    return is_leap_year(year) ? 29 : 28;
  default:
    return 0; // invalid month
  }
}

int days_from_start_of_year(int year, int month, int day) {
  int days = 0;
  for (int m = 1; m < month; ++m) {
    days += days_in_month(year, m);
  }
  days += day;
  return days;
}

int days_in_year(int year) {
  return is_leap_year(year) ? 366 : 365;
}

// Inclusive days between two dates: start and end (year, month, day)
int days_inclusive(int start_year, int start_month, int start_day,
                   int end_year, int end_month, int end_day) {
  if (start_year == end_year) {
    return days_from_start_of_year(end_year, end_month, end_day) - days_from_start_of_year(start_year, start_month, start_day) + 1;
  }

  int days = days_in_year(start_year) - days_from_start_of_year(start_year, start_month, start_day) + 1;
  for (int year = start_year + 1; year < end_year; ++year) {
    days += days_in_year(year);
  }
  days += days_from_start_of_year(end_year, end_month, end_day);
  return days;
}

//' @export
// [[Rcpp::export]]
Rcpp::List CASSIA_yearly(int start_year,
                         int end_year,

                         Rcpp::DataFrame weather,
                         std::vector<double> GPP_ref,

                         std::vector<double> pPREL,
                         Rcpp::DataFrame pCASSIA_parameters,
                         Rcpp::DataFrame pCASSIA_common,
                         Rcpp::DataFrame pCASSIA_ratios,
                         Rcpp::DataFrame pCASSIA_sperling,

                         std::vector<double> pPhydro,

                         int no_trees,

                         double needle_mass_in, // The value of this should be 0 if you want the needle value to be calculated
                         double Throughfall,

                         bool surplus_c,
                         bool nitrogen_change,
                         bool nitrogen_contrast,

                         double nitrogen_balance,

                         Rcpp::List settings) {

  /*
   * Read into structures
   */

  p1 parSite = make_p1(pPREL);
  p2 parGPP = make_p2(pPREL);
  p3 parET = make_p3(pPREL);
  p4 parSnowRain = make_p4(pPREL);
  p5 parWater = make_p5(pPREL);
  p7 parN = make_p7(pPREL);

  CASSIA_common common = make_common(pCASSIA_common);
  CASSIA_parameters parameters = make_CASSIA_parameters(pCASSIA_parameters, pCASSIA_sperling);
  CASSIA_ratios ratios = make_ratios(pCASSIA_ratios);

  Settings boolsettings = parseSettings(settings);
  phydro_canopy_parameters parPhydro = parPhydro_initalise(pPhydro);

  /*
   * Weather input made into vectors
   */

  weather_all climate = readWeatherVariables(weather, boolsettings.photosynthesis_as_input, boolsettings.preles, boolsettings.phydro);

  /*
   * Output structure set up
   */
  int simulation_time = days_inclusive(start_year, 1, 1, end_year, 12, 31);
  output_vector all_out{};
  initialize_output_vector(all_out, simulation_time);

  // TODO: add / change the original initialization!
  // Set up the initial conditions
  // days_per_year not defined as is later

  // Old
  double last_year_HH{275.4137};
  double last_year_maxN = 0.0;
  double GPP_mean = 0.0;
  double HH = 0.0;
  double GPP_sum_yesterday = 0.0;
  double GPP_sum = 0.0;
  double needles_last = 0.0;
  double winter_state = 1.0;
  std::vector<double> GPP_previous_sum(end_year - start_year + 2, parameters.GPP_initial);
  std::vector<double> potenital_growth_use(simulation_time, 0.0);
  std::vector<double> release(simulation_time, 0.0);

  /*
   * Structures set up
   */

  photosynthesis_out photosynthesis{};
  growth_state tree_state{};
  carbo_tracker sugar{};
  carbo_tracker starch{};
  carbo_tracker storage_term{};
  growth_out nitrogen_capacity{};
  uptake_structre uptake{};
  bool tree_alive = true;

  // Forward init values (previous day values) as first values of result vectors
  double CH = parameters.density_tree * parameters.carbon_share;
  double M_suc = 12 * common.M_C + 22 * common.M_H + 11 * common.M_O;
  double max_needles = 0.22;

  all_out.culm_growth.height[0] = parameters.h0;
  all_out.culm_growth.diameter[0] = parameters.D0;
  all_out.culm_growth.diameter_potential[0] = parameters.D0;
  all_out.culm_growth.roots[0] = 2.8; // Pauliina, 2019
  all_out.culm_growth.mycorrhiza[0] = 2.0; // TODO: make dynamic and find a better value

  repola_out repola_values{};
  if (needle_mass_in == 0) { // The value of this should be 0 if you want the needle value to be calculated
    repola_values = repola(all_out.culm_growth.diameter[0], all_out.culm_growth.height[0], parameters);
  } else {
    repola_values.needle_mass = needle_mass_in;
  }

  // Seed leaf_mass_correction's accumulating pool (actual_growth.cpp) at
  // the tree's initial Repola prediction, scaled by (n_age-1)/n_age - NOT
  // the full N_repo. Day 0 represents the state BEFORE this year's growth
  // has happened, same as any other day == 0: the fAPAR_and_LAI.cpp
  // formula itself reduces to exactly N_repo*(n_age-1)/n_age when
  // f_needle = 0 (no growth yet) - just the carryover of older retained
  // cohorts, before the newest one has grown in. Seeding at the full
  // N_repo overstated the tree's starting foliage. Mirrors height/
  // diameter's own explicit day-0 seeding (h0/D0 above) rather than
  // starting from an unseeded 0.0. Harmless when leaf_mass_correction =
  // FALSE: the old fAPAR_and_LAI.cpp pathway never reads index 0 (it is
  // memoryless).
  all_out.culm_growth.leaf_mass[0] = repola_values.needle_mass * (parameters.n_age - 1.0) / parameters.n_age;

  // Bootstraps C_N_year's and litter's first (day == 0, year 1)
  // computation in growth.cpp/actual_growth.cpp from this same seed,
  // rather than the struct default 0.0.
  tree_state.max_leaf_mass_running = all_out.culm_growth.leaf_mass[0];

  // Matching bootstrap for needle litter's own one-year-lagged normalizer
  // (growth.cpp's sum_degree_fL_last_year, function_structures.h), which
  // also defaults to 0.0 - without this, year 1 gets real growth (from
  // the fN bootstrap above) but zero litter at all (nothing to normalize
  // against), so the pool overshoots N_repo unopposed through years 1-2.
  // That overshoot then feeds makela_loss_boost_factor's tau_N=5 exponent
  // in year 3 - even a mild ~24% overshoot raised to the 5th power
  // produces a litter total ~2.9x a normal year's, crashing the pool to
  // near-zero and permanently locking it into a much-lower-than-N_repo
  // limit cycle it can never climb back out of (confirmed via debug
  // trace, 1998-2021: annual litter/target ratio is 0.00 (1998), 1.08
  // (1999), 2.86 (2000), then ~1.0 every year after). Unlike fN (a
  // pre-existing, independently-motivated phenology curve with its own
  // analytic-vs-empirical discussion above), sum_degree_fL_running's
  // annual total is remarkably stable and purely a property of the
  // fS^2/f_L/sLc mechanism itself (not weather-dependent the way fN's
  // g(T) sigmoid is) - traced directly across 24 simulated years
  // (1998-2021, debug trace): settles at 2.89-3.07, mean 2.99, every
  // single year. 3.0 used here as the same kind of bare empirical
  // constant as fN's 48.0, not independently calibrated.
  tree_state.sum_degree_fL_last_year = 3.0;

  // Organ-level sugar model uses per-organ initial values.
  // Non-organ-level-sugar model uses a single needle pool seeded from the aggregate sugar0/starch0.
  if (boolsettings.organ_level_sugar) {
    all_out.starch_vector.needles[0] = starch.needles = parameters.starch_needles00 = parameters.starch_needles0;
    all_out.sugar_vector.needles[0]  = sugar.needles  = parameters.sugar_needles00  = parameters.sugar_needles0;
    all_out.starch_vector.phloem[0]  = starch.phloem  = parameters.starch_phloem00  = parameters.starch_phloem0;
    all_out.sugar_vector.phloem[0]   = sugar.phloem   = parameters.sugar_phloem00   = parameters.sugar_phloem0;
    all_out.starch_vector.roots[0]   = starch.roots   = parameters.starch_roots00   = parameters.starch_roots0;
    all_out.sugar_vector.roots[0]    = sugar.roots    = parameters.sugar_roots00    = parameters.sugar_roots0;
    all_out.starch_vector.xylem_sh[0] = starch.xylem_sh = parameters.starch_xylem_sh00 = parameters.starch_xylem_sh0;
    all_out.sugar_vector.xylem_sh[0]  = sugar.xylem_sh  = parameters.sugar_xylem_sh00  = parameters.sugar_xylem_sh0;
    all_out.starch_vector.xylem_st[0] = starch.xylem_st = parameters.starch_xylem_st00 = parameters.starch_xylem_st0;
    all_out.sugar_vector.xylem_st[0]  = sugar.xylem_st  = parameters.sugar_xylem_st00  = parameters.sugar_xylem_st0;
  } else {
    // Single-pool model: total sugar and starch live in the needles slot
    all_out.sugar_vector.needles[0]  = sugar.needles  = parameters.sugar00  = parameters.sugar0;
    all_out.starch_vector.needles[0] = starch.needles = parameters.starch00 = parameters.starch0;
    // Pre-populate scalar totals for day 0 (sugar_model is skipped that day)
    all_out.sugar[0]  = sugar.needles;
    all_out.starch[0] = starch.needles;
  }

  // Nitrogen capacity and storage term start at 1 (non-limiting) before the
  // sugar model has run for the first time.
  nitrogen_capacity.needles = nitrogen_capacity.bud = nitrogen_capacity.wall =
    nitrogen_capacity.height = nitrogen_capacity.roots = 1.0;
  storage_term.needles = storage_term.phloem = storage_term.roots =
    storage_term.xylem_sh = storage_term.xylem_st = 1.0;

  // log_sugar() (and therefore storage_term_vector) is never reached on day 0
  // of the first year, since sugar_model() is skipped that day (see below) -
  // pre-populate it here so the "storage" output starts at 1, not 0.
  all_out.storage_term_vector.needles[0] = all_out.storage_term_vector.phloem[0] =
    all_out.storage_term_vector.roots[0] = all_out.storage_term_vector.xylem_sh[0] =
    all_out.storage_term_vector.xylem_st[0] = 1.0;

  all_out.nitrogen_capacity_vector.needles[0] = 1;
  all_out.nitrogen_capacity_vector.diameter[0] = 1;
  all_out.nitrogen_capacity_vector.height[0] = 1;
  all_out.nitrogen_capacity_vector.bud[0] = 1;
  all_out.nitrogen_capacity_vector.roots[0] = 1;

  // Needle litter's skew-normal xi (growth.cpp/actual_growth.cpp) needs
  // LAST year's fS<0.5 crossing day, which doesn't exist for the very
  // first simulated year. The weather driving fS is already fully known
  // in advance (not real look-ahead - this is the same complete,
  // historical input series the main loop below reads from, just read
  // once more here first), so start_year's own crossing day can be
  // found with a cheap, self-contained pre-pass (fS only ever depends on
  // TAir and its own running state - fS_model(), gpp.c - never on
  // growth/sugar/nitrogen state) and used to seed
  // tree_state.fS_cross_day_last_year before the main loop begins. That
  // makes the first simulated year behave exactly like every later one
  // from day 1 onward, with no bootstrap special-casing needed in
  // growth.cpp/actual_growth.cpp at all. Same -900 missing-data fixup
  // and fS_temp_sum=20.0 baseline the main loop itself uses (CASSIA.cpp
  // below), so this is deterministically identical to what the main
  // loop would compute for these same days anyway.
  if (boolsettings.leaf_mass_correction) {
    double fS_temp_sum_prepass = 20.0;
    int days_in_start_year = leap_year(start_year);
    for (int d = 0; d < days_in_start_year; d++) {
      if (boolsettings.preles && d > 0 && climate.TAir[d] < -900) {
        climate.TAir[d] = climate.TAir[d - 1];
      }
      double fS_prepass = fS_model(&fS_temp_sum_prepass, climate.TAir[d], parGPP);
      double doy_prepass = static_cast<double>(d) + 1.0;
      if (doy_prepass > 150 && fS_prepass < 0.5) {
        tree_state.fS_cross_day_last_year = doy_prepass;
        break;
      }
    }
  }

  /*
   * YEAR LOOP
   */
  int days_gone = 0;
  double LAI = 5; // TODO: check value

  for (int year = start_year; year <= end_year; year++) {

    /*
     * Daily output
     */

    // photosynthesis.fS itself is no longer reset here (2026-11-10) - it's
    // now computed fresh every day inside the day loop below, before the
    // photosynthesis dispatch, regardless of which model is active (see
    // fS_decoupled_from_preles_implementation_plan.md). The persistent
    // state behind it (tree_state.fS_temp_sum) carries its own reset
    // logic there, tied to soil_water_reset - so nothing here needs to
    // zero it first.
    bool fS_reached_one = false;
    bool fN_reached_one = false;

    /*
     * Yearly initialization
     */

    // B0, D00 and h00
    double B0 = M_PI/4.0 * pow(parameters.D0, 2.0);

    // TODO: should make this before the yearly thing?
    double D00 = parameters.D0; // Used to make the storage things in CASSIA
    double h00 = parameters.h0; // TODO: this has never been used
    if (boolsettings.xylogensis_option) {
      double LH0 = parameters.h_increment / (0.5 * parameters.sHc);
      double LN0 = parameters.n_length / (0.5 * parameters.sNc); // TODO: check the 2.64 which is in CASSIA_soil
      double LR0 = 2.0 * parameters.m_R_tot / parameters.sRc; // TODO: check if these parameters are updated somewhere
    }

    /*
     * NEEDLE MASS CALCULATION
     *
     * If it is the first year, then there is no needle mass initialized,
     * after this if there is growth in the model the needle mass is calculated based on last year
     * if there is no growth then the needle mass stays at the originally calculated value
     */

    /*
     * Indexes
     */

    int index = days_gone - 1;
    if (index < 0) {
      index = 0;
    }

    if (boolsettings.needle_mass_grows) {
      repola_values = repola(all_out.culm_growth.diameter[index],
                             all_out.culm_growth.height[index], parameters); // Needle mass is then calculated on the net D and h values
    } else {
      repola_values.needle_mass = needle_mass_in;
    }

    /*
     * Days per year
     */

    int days_per_year = leap_year(year);

    /*
     * DAYS LOOP
     */
    for (int day = 0; day < days_per_year; day++) {
      all_out.year[days_gone + day] = year;
      all_out.day[days_gone + day] = day + 1;

      /*
       * Weather checks!
       */

      if (boolsettings.preles) {
        if (day > 0) {
          if (climate.PAR[day + days_gone] < -900) climate.PAR[day + days_gone] = climate.PAR[day + days_gone -1];
          if (climate.TAir[day + days_gone] < -900) climate.TAir[day + days_gone] = climate.TAir[day + days_gone - 1];
          if (climate.VPD[day + days_gone] < 0 || climate.VPD[day + days_gone] > 6) climate.VPD[day + days_gone] = climate.VPD[day + days_gone -1];
          if (climate.Precip[day + days_gone] < 0) climate.Precip[day + days_gone] = climate.Precip[day + days_gone -1] * 0.3;
          /* On avg. P+1=0.315*P
           * (in Sodis & Hyde) */
          if (climate.CO2[day + days_gone] < 0) climate.CO2[day + days_gone] = climate.CO2[day + days_gone -1];
        }
      }

      /*
       * LAI
       *
       * There are yearly, but not daily dependencies other than environmental states here!
       */

      // At the start, have a think about this! Isn't it also in tree state?
      compute_fAPAR_used(day,
                         days_gone,
                         LAI,
                         max_needles,
                         fS_reached_one,
                         fN_reached_one,
                         repola_values,
                         boolsettings,
                         parameters,
                         tree_state,
                         photosynthesis,
                         climate,
                         all_out);

      /*
       * Photosynthesis
       */

      // fS (Mäkelä et al. 2004 seasonality model): computed here,
      // unconditionally, once per day, for every photosynthesis pathway -
      // not just preles. reset_state: soil_water_reset=TRUE (default)
      // fires this every year, identical to the old unconditional day==0
      // behavior preles_cpp() used to apply only to itself; FALSE fires it
      // only on the true first day of the whole run (year==start_year),
      // letting it carry over continuously across years instead - needed
      // for ORCHIMIC's continuous multi-year runs
      // (manscript_orchemic/PROJECT_INSTRUCTIONS.md). 20.0 baseline
      // matches PRELES's own Initials_snow.S (preles_test.cpp).
      bool preles_reset_state = (day == 0) && (boolsettings.soil_water_reset || year == start_year);
      if (preles_reset_state) {
        tree_state.fS_temp_sum = 20.0;
      }
      photosynthesis.fS = fS_model(&tree_state.fS_temp_sum, climate.TAir[days_gone + day], parGPP);

      double photosynthesis_per_stem = 0.0;
      if (boolsettings.photosynthesis_as_input) {
        photosynthesis.GPP = climate.Photosynthesis_IN[day + days_gone];
        photosynthesis.ET = 0.0;
        photosynthesis.SoilWater = 0.0;
        photosynthesis_per_stem = climate.Photosynthesis_IN[day + days_gone] / no_trees * 10000/1000;
      } else if (boolsettings.preles) {
        photosynthesis = preles_cpp(day, preles_reset_state, climate.PAR[days_gone + day], climate.TAir[days_gone + day], climate.Precip[days_gone + day],
                                    climate.VPD[days_gone + day], climate.CO2[days_gone + day], all_out.fAPAR[days_gone + day], photosynthesis.fS,
                                    parSite, parGPP, parET, parSnowRain, parWater, 0.0, 1);
        photosynthesis_per_stem = photosynthesis.GPP / no_trees * 10000/1000;
      } else if (boolsettings.phydro) {
        // Wired in 2026-10-29: calc_plant_assimilation_rate() and parPhydro_initalise()
        // already existed (used by the disabled CASSIA_eeo export) but this branch was
        // simply missing from the main dispatch, which is why phydro=TRUE silently gave
        // GPP=0 before. crown_area formula and the zeta proxy are adapted from the
        // commented-out ecoevolutionary.cpp reference (which itself flags "TODO:
        // something wrong here!" on the crown_area line) - first working pass, not
        // verified against real canopy geometry.
        int index_ref = days_gone + day - 1;
        if (index_ref < 0) index_ref = 0;

        double height_prev = all_out.culm_growth.height[index_ref];
        double diameter_prev = all_out.culm_growth.diameter[index_ref];
        double crown_area = M_PI * 6000 / (4 * 75) * height_prev * diameter_prev;

        // Phydro-local LAI, kept separate from the shared `LAI` (=5, used by
        // compute_fAPAR_used()/preles above with its own 0.52-coefficient
        // fAPAR formula - changing that shared constant would alter the
        // preles/photosynthesis_as_input paths too). 4 is the real,
        // site-calibrated value for Hyytiala ("lai0 = 4 # (Hyytiala data)" in
        // the Plant-FATE/TreeLife companion codebase's boreal config),
        // confirmed 2026-10-29. calc_plant_assimilation_rate() used to ignore
        // whatever lai it was given anyway (see the fix to PlantFateLogic.cpp
        // removing its internal `lai = 1.8` Amazon-borrowed override).
        double LAI_phydro = 4.0;
        double zeta = LAI_phydro / all_out.culm_growth.roots[index_ref];

        parPhydro.tau_weather = 7;
        parPhydro.dt = (day == 0) ? 0.0 : 1.0 / days_per_year;

        // Unit conversions confirmed 2026-10-29 against the Plant-FATE/TreeLife
        // companion codebase (~/Documents/Austria/Plant-FATE/src/climate_stream.cpp),
        // which applies the same two conversions to the same raw ERA5-style input
        // before calling phydro: VPD is stored in hPa here (ERA5_CASSIA_phydro.csv
        // ranges 0-20.6, too high for kPa) and must become Pa; SWP is stored as a
        // positive magnitude ("-MPa" with the sign already stripped, range
        // 0.017-0.133) and must be negated to a true soil water potential. Without
        // these the phydro solver failed to converge on 364/365 days in testing and
        // the one "successful" day produced a physically impossible GPP (~-2.5e23).
        double VPD_pa = climate.VPD[days_gone + day] * 100.0;   // hPa -> Pa
        double SWP_mpa = climate.SWP[days_gone + day] * (-1.0); // positive magnitude -> true (negative) MPa

        // The phydro solver (LBFGS line search) occasionally fails to converge,
        // most often in cold weather. Rather than aborting the whole run, fall
        // back to the previous day's photosynthesis values and keep going -
        // photosynthesis.GPP/ET already hold the last successfully-computed
        // day's values at this point, so simply skipping the overwrite below
        // achieves the carry-forward. A warning with the day/year is printed
        // so the frequency of these misfires can be counted afterward.
        try {
          PlantAssimilationResult photosynthesis_phydro = calc_plant_assimilation_rate(
            climate.PAR[days_gone + day], climate.PAR_max[days_gone + day], climate.TAir[days_gone + day],
            VPD_pa, climate.Precip[days_gone + day], climate.CO2[days_gone + day],
            climate.Nitrogen[days_gone + day], climate.PA[days_gone + day], SWP_mpa,
            parPhydro, LAI_phydro, crown_area, height_prev, zeta, day);

          // The solver can also return "successfully" (no exception) with a
          // NaN/Inf result - e.g. a 0/0 or similar inside the optimiser that
          // doesn't itself throw. Found 2026-10-29: this happened silently on
          // day 4 of a 2015 test run, and because photosynthesis.GPP feeds a
          // cumulative smoothing term (state.S_GPP in growth.cpp) that is
          // never reset, a single silent NaN here permanently poisons every
          // later day, surfacing ~77 days later as a diameter/crown_area NaN
          // crash with no clue it started here. Treat non-finite output the
          // same as a thrown exception: warn and keep the previous day's
          // values rather than letting it propagate.
          if (!std::isfinite(photosynthesis_phydro.gpp) || !std::isfinite(photosynthesis_phydro.trans)) {
            Rcpp::Rcout << "Warning: phydro solver returned a non-finite result on day " << (day + 1)
                        << " of year " << year << " (gpp=" << photosynthesis_phydro.gpp
                        << ", trans=" << photosynthesis_phydro.trans
                        << ") - using previous day's photosynthesis values.\n";
          } else {
            photosynthesis.GPP = photosynthesis_phydro.gpp;
            photosynthesis.ET = photosynthesis_phydro.trans;
            photosynthesis.SoilWater = 0.0;
          }
        } catch (const std::exception& e) {
          Rcpp::Rcout << "Warning: phydro solver did not converge on day " << (day + 1)
                      << " of year " << year << " (" << e.what()
                      << ") - using previous day's photosynthesis values.\n";
        }
        photosynthesis_per_stem = photosynthesis.GPP / no_trees * 10000/1000;
      } else {
        std::cout << "There is no photosynthesis model chosen!\n";
      }

      log_photosynthesis(day, days_gone, photosynthesis, all_out);

      /*
       * Potential Growth
       *
       * In terms of the adaptation from the R code, the potential values are not altered by daily processes so still calculate them for a year
       */

      double day_gpp = day;
      if (day_gpp > 364) day_gpp = 364;
      growth(day, days_gone, year,
             tree_state,
             all_out,
             climate.TAir[days_gone + day],
             climate.TSoil_A[days_gone + day],
             climate.TSoil_B[days_gone + day],
             climate.Soil_Moisture[days_gone + day],
             photosynthesis.GPP,
             GPP_ref[day_gpp],
             boolsettings,
             common, parameters, ratios,
             CH,
             B0,
             GPP_mean,
             GPP_previous_sum[year-start_year],
             // Last iteration value
             last_year_HH,
             days_per_year,
             repola_values.needle_mass);

      // Saved for the next iteration
      release[days_gone + day] = tree_state.use;
      int lim = static_cast<int>(std::ceil(parameters.tau_Ee));
      if (day > (lim - 1)) {
        tree_state.release = release[days_gone + day - lim];
      } else {
        tree_state.release = 0.0;
      }

      /*
       * Respiration
       *
       * There are yearly, but not daily dependencies other than weather conditions here!
       */

      respiration(day, tree_state, all_out, parameters, ratios,
                  climate.TAir[days_gone + day], climate.TSoil_A[days_gone + day],
                  boolsettings,
                  B0);

      /*
       * Sugar
       */

      if (!(day == 0 && year == start_year)) {
        // Note: mycorrhiza is always passive at the moment
        // TODO: need to add some more boolsettings!
        sugar_model(day,
                    days_gone,
                    climate.TAir[days_gone + day],
                    climate.PAR[days_gone + day],
                    photosynthesis_per_stem,
                    parameters,
                    common,
                    nitrogen_change,
                    nitrogen_contrast,
                    TRUE,
                    surplus_c,
                    tree_alive,
                    boolsettings.organ_level_sugar,
                    boolsettings.storage_reset,
                    tree_state,
                    nitrogen_balance,
                    winter_state,
                    uptake,
                    sugar,
                    starch,
                    storage_term,
                    nitrogen_capacity,
                    all_out);
      }

      /*
       * Actual growth
       */

      actual_growth(day, days_gone,
                    parameters, common,
                    storage_term, photosynthesis,
                    tree_state,
                    all_out,
                    boolsettings,
                    nitrogen_capacity);

      /*
       * Updating some parameters
       */

      // TODO: what does this do?
      HH = tree_state.HH;
      potenital_growth_use[days_gone + day] = tree_state.use;

      GPP_sum_yesterday = GPP_sum;

      //if (day == 0 || day == days_per_year-1) {
        //all_out.sugar_vector.print(day + days_gone, "sugar");
        //all_out.starch_vector.print(day + days_gone, "starch");
        //all_out.culm_growth.print(day + days_gone, "culm");
      //}

    } // End of the days loop

    // Updating parameters!

    last_year_HH = HH;

    GPP_previous_sum[year - start_year + 1] = std::accumulate(all_out.photosynthesis.GPP.begin() + days_gone + 182,
                                                              all_out.photosynthesis.GPP.begin() + days_gone + 245 + 1, 0.0);

    double start_needles = all_out.culm_growth.needles[days_gone];
    double end_needles   = all_out.culm_growth.needles[days_gone + days_per_year - 1];
    max_needles          = end_needles - start_needles;

    days_gone = days_gone + days_per_year;
  } // End of the years loop


  ///////////////////////
  // Output
  ///////////////////////

  // Dataframes

  // printOutputVectorSizes(all_out);

  if (nitrogen_change) {
    Rcpp::DataFrame df5 = Rcpp::DataFrame::create(
      Rcpp::_["root_uptake"] = Rcpp::wrap(all_out.uptake_vector.root_uptake),
      Rcpp::_["ectomycorrhizal_uptake"] = Rcpp::wrap(all_out.uptake_vector.ectomycorrhizal_uptake),
      Rcpp::_["total_uptake"] = Rcpp::wrap(all_out.uptake_vector.total_uptake),
      Rcpp::_["ectomycorrhizal_transfer"] = Rcpp::wrap(all_out.uptake_vector.ectomycorrhizal_transfer)
    );

    return Rcpp::List::create(
      Rcpp::_["Growth"] = createGrowthDataFrame(all_out),
      Rcpp::_["Sugar"] = createSugarDataFrame(all_out),
      Rcpp::_["Preles"] = createPrelesDataFrame(all_out),
      Rcpp::_["Culm_Growth"] = createCulmGrowthDataFrame(all_out),
      Rcpp::_["Turnover"] = createTurnoverDataFrame(all_out),
      Rcpp::_["Uptake"] = df5
    );

  } else {
    return Rcpp::List::create(
      Rcpp::_["Growth"] = createGrowthDataFrame(all_out),
      Rcpp::_["Sugar"] = createSugarDataFrame(all_out),
      Rcpp::_["Preles"] = createPrelesDataFrame(all_out),
      Rcpp::_["Culm_Growth"] = createCulmGrowthDataFrame(all_out),
      Rcpp::_["Turnover"] = createTurnoverDataFrame(all_out)
    );

  }

}
