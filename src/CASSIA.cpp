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

  /*
   * Weather input made into vectors
   */

  weather_all climate = readWeatherVariables(weather, boolsettings.photosynthesis_as_input, boolsettings.preles, boolsettings.phydro);

  /*
   * Output structure set up
   */
  int simulation_time = days_inclusive(start_year, 1, 1, end_year, 12, 31); // TODO: create a function that creates this size output
  output_vector all_out{};
  initialize_output_vector(all_out, simulation_time);


  // TODO: add / change the original initialization!
  // Set up the initial conditions
  // days_per_year not defined as is later

  // Old
  needle_cohorts last_cohorts{};
  double last_year_HH{};
  double last_year_maxN{};
  double GPP_mean{};
  std::vector<double> GPP_previous_sum = {parameters.GPP_initial};
  std::vector<double> potenital_growth_use{};

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
  double max_needles = 0.6;

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

  all_out.starch_vector.needles[0] = starch.needles = parameters.starch_needles00 = parameters.starch_needles0;
  all_out.sugar_vector.needles[0] = sugar.needles = parameters.sugar_needles00 = parameters.sugar_needles0;
  all_out.starch_vector.phloem[0] = starch.phloem = parameters.starch_phloem00 = parameters.starch_phloem0;
  all_out.sugar_vector.phloem[0] = sugar.phloem = parameters.sugar_phloem00 = parameters.sugar_phloem0;
  all_out.starch_vector.roots[0] = starch.roots = parameters.starch_roots00 = parameters.starch_roots0;
  all_out.sugar_vector.roots[0] = sugar.roots = parameters.sugar_roots00 = parameters.sugar_roots0;
  all_out.starch_vector.xylem_sh[0] = starch.xylem_sh = parameters.starch_xylem_sh00 = parameters.starch_xylem_sh0;
  all_out.sugar_vector.xylem_sh[0] = sugar.xylem_sh = parameters.sugar_xylem_sh00 = parameters.sugar_xylem_sh0;
  all_out.starch_vector.xylem_st[0] = starch.xylem_st = parameters.starch_xylem_st00 = parameters.starch_xylem_st0;
  all_out.sugar_vector.xylem_st[0] = sugar.xylem_st = parameters.sugar_xylem_st00 = parameters.sugar_xylem_st0;

  all_out.nitrogen_capacity_vector.needles[0] = 1;
  all_out.nitrogen_capacity_vector.diameter[0] = 1;
  all_out.nitrogen_capacity_vector.height[0] = 1;
  all_out.nitrogen_capacity_vector.bud[0] = 1;
  all_out.nitrogen_capacity_vector.roots[0] = 1;

  /*
   * YEAR LOOP
   */
  int days_gone = 0;
  double LAI = 5; // TODO: check value

  for (int year = start_year; year <= end_year; year++) {
    /*
     * Daily output
     */

    // TODO: think about this bit!
    photosynthesis.fS = 0.0;
    bool fS_reached_one = false;
    bool fN_reached_one = false;
    std::vector<double> release;

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
      repola_values = repola(all_out.culm_growth.diameter[index], all_out.culm_growth.height[index], parameters); // Needle mass is then calculated on the next D0 and h0 values
    } else {
      repola_values.needle_mass = needle_mass_in;
    }

    // NOTE: not used anywhere!
    needle_cohorts needles_cohorts;
    if (year > start_year) {
      needles_cohorts.year_1 = 31.24535 / parameters.n_length * repola_values.needle_mass / 3.0;
      needles_cohorts.year_2 = last_cohorts.year_1;
      needles_cohorts.year_3 = last_cohorts.year_2;
    }

    double HH{}, GPP_sum_yesterday{}, GPP_sum{};
    double needles_last{};
    if (year == start_year) {
      last_year_HH = 275.4137;
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
          if (climate.TAir[day + days_gone] < -900) climate.TAir[day + days_gone] = climate.TAir[day + days_gone- 1];
          if (climate.VPD[day + days_gone] < 0 || climate.VPD[day + days_gone] > 6) climate.VPD[day + days_gone] = climate.VPD[day + days_gone -1];
          if (climate.Precip[day + days_gone] <    0) climate.Precip[day + days_gone] = climate.Precip[day + days_gone -1] * 0.3;
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

      double photosynthesis_per_stem = 0.0;
      if (boolsettings.photosynthesis_as_input) {
        photosynthesis.GPP = climate.Photosynthesis_IN[day + days_gone];
        photosynthesis.ET = 0.0;
        photosynthesis.SoilWater = 0.0;
        photosynthesis_per_stem = climate.Photosynthesis_IN[day + days_gone] / 1010 * 10000/1000;
      } else if (boolsettings.preles) {
        photosynthesis = preles_cpp(day, climate.PAR[days_gone + day], climate.TAir[days_gone + day], climate.Precip[days_gone + day],
                                    climate.VPD[days_gone + day], climate.CO2[days_gone + day], all_out.fAPAR[days_gone + day],
                                    parSite, parGPP, parET, parSnowRain, parWater, 0.0, 1);
        photosynthesis_per_stem = photosynthesis.GPP / 1010 * 10000/1000;
      } else {
        std::cout << "There is no photosynthesis model chosen!\n";
      }

      log_photosynthesis(day, days_gone, photosynthesis, all_out);

      /*
       * Potential Growth
       *
       * In terms of the adaptation from the R code, the potential values are not altered by daily processes so still calculate them for a year
       */

      growth(day, days_gone, year,
             tree_state,
             all_out,
             climate.TAir[days_gone + day],
             climate.TSoil_A[days_gone + day],
             climate.TSoil_B[days_gone + day],
             climate.Soil_Moisture[days_gone + day],
             photosynthesis.GPP,
             GPP_ref[days_gone + day],
             boolsettings,
             common, parameters, ratios,
             CH,
             B0,
             GPP_mean,
             GPP_previous_sum[year-start_year],
             // Last iteration value
             last_year_HH,
             days_per_year);

      // Saved for the next iteration
      release.push_back(tree_state.en_pot_growth);
      double lim = std::ceil(parameters.tau_Ee);
      if (day > (lim-1)) {
        tree_state.release = release[day-lim];
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
                    boolsettings.sperling_model,
                    tree_state,
                    nitrogen_balance,
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
      potenital_growth_use.push_back(tree_state.use);

      GPP_sum_yesterday = GPP_sum;
    } // End of the days loop

    // Updating parameters!

    last_year_HH = HH;

    GPP_previous_sum.push_back(std::accumulate(all_out.photosynthesis.GPP.begin() + days_gone + 182, all_out.photosynthesis.GPP.begin() + days_gone + 245 + 1, 0.0));

    double start_needles = all_out.culm_growth.needles[days_gone];
    double end_needles   = all_out.culm_growth.needles[days_gone + days_per_year - 1];
    double max_needles   = end_needles - start_needles;

    std::cout << " max_needles " << max_needles << "\n";

    days_gone = days_gone + days_per_year;

    last_cohorts.year_1 = needles_cohorts.year_1; // TODO: currently the growth doesn't really have an effect on this - should it?
    last_cohorts.year_2 = needles_cohorts.year_2;
    last_cohorts.year_3 = needles_cohorts.year_3;
    // TODO: make it possible to make this for 5 years?
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
      Rcpp::_["Uptake"] = df5
    );

  } else {
    return Rcpp::List::create(
      Rcpp::_["Growth"] = createGrowthDataFrame(all_out),
      Rcpp::_["Sugar"] = createSugarDataFrame(all_out),
      Rcpp::_["Preles"] = createPrelesDataFrame(all_out),
      Rcpp::_["Culm_Growth"] = createCulmGrowthDataFrame(all_out)
    );

  }

}
