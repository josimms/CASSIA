#include "CASSIA.h"

int leap_year(int year) {
  if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) {
    return 366;
  } else {
    return 365;
  }
};

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
   * Structures set up
   */

  // Forward init values (previous day values) as first values of result vectors

  bool tree_alive = true;

  double CH = parameters.density_tree * parameters.carbon_share;
  double M_suc = 12 * common.M_C + 22 * common.M_H + 11 * common.M_O;
  double max_needles = 0.03;

  carbo_tracker Ad;
  Ad.needles = parameters.Ad0_needles;
  Ad.phloem = parameters.Ad0_phloem;
  Ad.xylem_sh = parameters.Ad0_xylem_sh;
  Ad.xylem_st = parameters.Ad0_xylem_st;
  Ad.roots = parameters.Ad0_roots;

  parameters.starch_needles00 = parameters.starch_needles0;
  parameters.sugar_needles00 = parameters.sugar_needles0;
  parameters.starch_phloem00 = parameters.starch_phloem0;
  parameters.sugar_phloem00 = parameters.sugar_phloem0;
  parameters.starch_roots00 = parameters.starch_roots0;
  parameters.sugar_roots00 = parameters.sugar_roots0;
  parameters.starch_xylem_sh00 = parameters.starch_xylem_sh0;
  parameters.sugar_xylem_sh00 = parameters.sugar_xylem_sh0;
  parameters.starch_xylem_st00 = parameters.starch_xylem_st0;
  parameters.sugar_xylem_st00 = parameters.sugar_xylem_st0;


  repola_out repola_values;
  if (needle_mass_in == 0) { // The value of this should be 0 if you want the needle value to be calculated
    repola_values = repola(parameters); // TODO: fix
  } else {
    repola_values.needle_mass = needle_mass_in;
  }

  /*
   * Vectors between iterations
   */

  growth_values_out growth_values_for_next_iteration = growth_values_out_init();
  carbo_balance sugar_values_for_next_iteration = carbo_balance_init();
  ring_width_out previous_ring_width = ring_width_out_init();

  double height_next_year = parameters.h0;
  double diameter_next_year = 1000 * parameters.D0;
  double diameter_potential_next_year = 1000 * parameters.D0;
  double roots_next_year = 2.0; // TODO: make dynamic
  double mycorrhiza_next_year = 2.0; // TODO: make dynamic
  double xylem_edge_estimate, xylem_edge_estimate_internal;

  /*
   * Vectors for the outputs
   */

  growth_vector potential_growth_output;
  growth_vector actual_growth_output;
  sugar_values_vector sugar_values_output;
  resp_vector respiration_output;
  growth_vector culm_growth;
  growth_vector culm_growth_internal;
  needle_cohorts last_cohorts;
  double last_year_HH;
  double last_year_maxN;
  double GPP_mean;
  std::vector<double> GPP_previous_sum;
  GPP_previous_sum.push_back(parameters.GPP_initial);
  double respiration_maintanence;
  std::vector<double> potenital_growth_use;

  photosynthesis_out photosynthesis;
  photo_out_vector photosynthesis_output;

  uptake_structre_vector uptake_values_output;

  /*
   * YEAR LOOP
   */
  std::vector<int> years;
  std::vector<int> days;

  std::vector<int> years_for_runs;
  int years_temparary = start_year;
  double end = 2*(end_year - start_year + 1);
  for (int i=1; i <= end; i++) {
    years_for_runs.push_back(years_temparary);
    if (i%2!=0) {
      years_temparary = years_temparary;
    } else {
      years_temparary = years_temparary+1;
    }
  }

  int final_year = 1;
  int days_gone = 0;

  double LAI = 5;
  std::vector<double> LAI_within_year_vector;

  // Temperature equilibrium for the sugar model
  //	# Compute initial Te by the mean temperature for the first week of # Semptemver plus 3C (for the exponential nature of the curves), original was October
  carbo_tracker equilibrium_temperature = carbo_tracker_init();
  double equilibrium_temperature_init = (climate.TAir[244] + climate.TAir[245] + climate.TAir[246] + climate.TAir[247] + climate.TAir[248] + climate.TAir[249] + climate.TAir[250] + climate.TAir[251]) / 7 + 3;
  equilibrium_temperature.needles = equilibrium_temperature.phloem = equilibrium_temperature.roots = equilibrium_temperature.xylem_sh = equilibrium_temperature.xylem_st = equilibrium_temperature_init;

  for (int year : years_for_runs)  {
    /*
     * Daily output
     */

    carbo_tracker carbo_tracker_vector;
    xylogensis_out xylogensis_vector;
    photosynthesis.fS = 0.0;

    std::vector<double> release;

    /*
     * Yearly initialization
     */

    // B0, D00 and h00
    double B0 = M_PI/4.0 * pow(parameters.D0, 2.0); // TODO: B0 is the start of Bud burst right?
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

    if (boolsettings.needle_mass_grows) {
      repola_values = repola(parameters); // Needle mass is then calculated on the next D0 and h0 values
    } else {
      repola_values.needle_mass = needle_mass_in;
    }

    needle_cohorts needles_cohorts;
    if (year > start_year) {
      // TODO: parameters, add the year 5
      needles_cohorts.year_1 = 31.24535 / parameters.n_length * repola_values.needle_mass / 3.0;
      needles_cohorts.year_2 = last_cohorts.year_1;
      needles_cohorts.year_3 = last_cohorts.year_2;
    }
    double HH, GPP_sum_yesterday, GPP_sum;;
    double needles_last;
    if (year == start_year) {
      last_year_HH = 275.4137;
    }

    // TODO: There should also be xylogenesis dependency here!

    /*
     * Days per year
     */

    int days_per_year = leap_year(year);
    /*
     * Yearly initial conditions updated
     */

    // Set up the initial conditions
    yearly_in yearly = yearly_initial_conditions(days_per_year); // TODO: growth should be added to this!

    /*
     * DAYS LOOP
     */
    int weather_index;
    photosynthesis_out photosynthesis_old;
    for (int day = 0; day < days_per_year; day++) {
      // std::cout << "Year " << year << " Day " << day;

      weather_index = days_gone + day;

      /*
       * Weather checks!
       */

      if (boolsettings.preles) {
        if (weather_index > 0) {
          if (climate.PAR[weather_index] < -900) climate.PAR[weather_index] = climate.PAR[weather_index-1];
          if (climate.TAir[weather_index] < -900) climate.TAir[weather_index] = climate.TAir[weather_index-1];
          if (climate.VPD[weather_index] < 0 || climate.VPD[weather_index] > 6) climate.VPD[weather_index] = climate.VPD[weather_index-1];
          if (climate.Precip[weather_index] <    0) climate.Precip[weather_index] = climate.Precip[weather_index-1] * 0.3;
          /* On avg. P+1=0.315*P
           * (in Sodis & Hyde) */
          if (climate.CO2[weather_index] < 0) climate.CO2[weather_index] = climate.CO2[weather_index-1];
        }
      }

      /*
       * PHOTOSYNTHESIS
       *
       * There are yearly, but not daily dependencies other than environmental states here!
       */

      double fAPAR_used = 0.0;
      double fS_out = 0.0;
      double needle_growth = 0.0;
      double LAI_within_year = 0.0;
      if (day > 0) {
        fS_out = photosynthesis_output.fS[weather_index-1];
        needle_growth = actual_growth_output.needles[weather_index-1];
      }
      if (!boolsettings.photosynthesis_as_input & boolsettings.fAPAR_Tian) {
        // Uses the method from Tian 2021
        // Extinction coefficient 0.52 is from Tian 2021 as well
        // LAI value is fairly constant if we look at Rautiainen 2012, LAI ~ 3
        double f_modifer = 0.0;
        if (max_needles == 0.0) {
          f_modifer = 0.0;
        } else {
          f_modifer = needle_growth/max_needles; // the actual growth divided by the maximum per year
        }
        if (day < 182) { // TODO: I decided that the start of July is the end of spring
          LAI_within_year = LAI*(parameters.n_age - 1)/parameters.n_age + (1/parameters.n_age)*f_modifer*LAI; // TODO: this varies too much!
        } else if (day > 244) { // TODO: I decided that the end of august is the start of autumn
          LAI_within_year = LAI*(parameters.n_age - 1)/parameters.n_age + (1/parameters.n_age)*fS_out*LAI;
        } else {
          LAI_within_year = LAI;
        }
        fAPAR_used = (1.0 - std::exp(-0.52 * LAI_within_year));  // TODO: Check this is sensible
        if (std::isnan(LAI_within_year)) {
          std::cout << "Day " << day << " LAI_within_year was NaN so replaced. TODO: fix\n";
          fAPAR_used = 0.7;
        }
        if (std::isnan(fAPAR_used)) {
          std::cout << "Day " << day << " fAPAR_used was NaN so replaced. TODO: fix\n";
          fAPAR_used = 0.7;
        }
      } else if (!boolsettings.photosynthesis_as_input & !boolsettings.fAPAR_Tian & boolsettings.preles) {
        fAPAR_used = climate.fAPAR[weather_index];
      } else {
        fAPAR_used = 0.0;
      }
      if (final_year%2!=0) {
        LAI_within_year_vector.push_back(LAI_within_year);
      }

      double photosynthesis_per_stem = 0.0;
      double GPP = 0.0;
      double ET = 0.0;
      double SoilWater = 0.0;
      if (boolsettings.photosynthesis_as_input) {
        GPP = photosynthesis.GPP = climate.Photosynthesis_IN[weather_index];
        ET = photosynthesis.ET = 0.0;
        SoilWater = photosynthesis.SoilWater = 0.0;
        photosynthesis_per_stem = climate.Photosynthesis_IN[weather_index] / 1010 * 10000/1000;
        if (final_year%2!=0) {
          photosynthesis_output.GPP.push_back(photosynthesis.GPP);
          photosynthesis_output.ET.push_back(photosynthesis.ET);
          photosynthesis_output.SoilWater.push_back(photosynthesis.SoilWater);
          photosynthesis_output.fAPAR.push_back(fAPAR_used);
        }
      } else if (boolsettings.preles) {
        if (final_year%2!=0) {
          photosynthesis = preles_cpp(weather_index, climate.PAR[weather_index], climate.TAir[weather_index], climate.Precip[weather_index],
                                      climate.VPD[weather_index], climate.CO2[weather_index], fAPAR_used,
                                      parSite, parGPP, parET, parSnowRain, parWater, 0.0, 1);
          photosynthesis_output.GPP.push_back(photosynthesis.GPP);
          photosynthesis_output.ET.push_back(photosynthesis.ET);
          photosynthesis_output.SoilWater.push_back(photosynthesis.SoilWater);
          photosynthesis_output.fS.push_back(photosynthesis.fS);
          photosynthesis_output.fAPAR.push_back(fAPAR_used);
          GPP = photosynthesis_output.GPP[weather_index];
          photosynthesis_per_stem = GPP / 1010 * 10000/1000;
          ET = photosynthesis_output.ET[weather_index];
          SoilWater = photosynthesis_output.SoilWater[weather_index];
        } else {
          GPP = photosynthesis_output.GPP[weather_index];
          photosynthesis_per_stem = GPP / 1010 * 10000/1000;
          ET = photosynthesis_output.ET[weather_index];
          SoilWater = photosynthesis_output.SoilWater[weather_index];
        }
      } else {
        std::cout << "There is no photosynthesis model chosen!\n";
      }

      /*
       * Potential Growth
       *
       * In terms of the adaptation from the R code, the potential values are not altered by daily processes so still calculate them for a year
       */

      growth_out potential_growth = growth(day, year, climate.TAir[weather_index], climate.TSoil_A[weather_index], climate.TSoil_B[weather_index], climate.Soil_Moisture[weather_index], GPP, GPP_ref[day],
                                           boolsettings.root_as_Ding, boolsettings.xylogensis_option, boolsettings.environmental_effect_xylogenesis, boolsettings.sD_estim_T_count,
                                           common, parameters, ratios,
                                           CH, B0, GPP_mean, GPP_previous_sum[year-start_year],
                                           boolsettings.LH_estim, boolsettings.LN_estim, boolsettings.LD_estim,
                                           boolsettings.tests,
                                           // Last iteration value
                                           growth_values_for_next_iteration, last_year_HH,
                                           days_per_year);
      if (final_year%2!=0) {
        potential_growth_output.height.push_back(potential_growth.height);
        potential_growth_output.needles.push_back(potential_growth.needles);
        potential_growth_output.roots.push_back(potential_growth.roots);
        potential_growth_output.diameter.push_back(potential_growth.diameter);
        potential_growth_output.bud.push_back(potential_growth.bud);
        potential_growth_output.g.push_back(potential_growth.g);
        potential_growth_output.en_pot_growth.push_back(potential_growth.use);
        potential_growth_output.ring_width.push_back(potential_growth.previous_values.pot_mm);
      }
      // Saved for the next iteration
      growth_values_for_next_iteration = potential_growth.previous_values;
      release.push_back(potential_growth.previous_values.en_pot_growth);
      double lim = std::ceil(parameters.tau_Ee);
      if (day > (lim-1)) {
        potential_growth.release = release[day-lim];
      } else {
        potential_growth.release = 0.0;
      }

      /*
       * Respiration
       *
       * There are yearly, but not daily dependencies other than weather conditions here!
       */

      respiration_out resp = respiration(day, parameters, ratios, repola_values,
                                         climate.TAir[weather_index], climate.TSoil_A[weather_index],
                                         boolsettings.temp_rise, boolsettings.Rm_acclimation, boolsettings.mN_varies,
                                         // parameters that I am not sure about
                                         B0);

      /*
       * Sugar
       */

      double root_mass = 0;
      double mycorrhizal_biomass = 0;
      double xylem_st_mass;
      double phloem_mass;
      if (day == 0) {
        root_mass = roots_next_year;
        mycorrhizal_biomass = mycorrhiza_next_year;

      } else {
        if (final_year%2!=0) {
          root_mass = culm_growth_internal.roots[weather_index-1];
          mycorrhizal_biomass = culm_growth_internal.mycorrhiza[weather_index-1];
          xylem_st_mass = culm_growth_internal.xylem_st[weather_index-1];
          phloem_mass = culm_growth_internal.phloem[weather_index-1];
        } else {
          root_mass = culm_growth.roots[weather_index-1];
          mycorrhizal_biomass = culm_growth.mycorrhiza[weather_index-1];
          xylem_st_mass = culm_growth.xylem_st[weather_index-1];
          phloem_mass = culm_growth.phloem[weather_index-1];
        }
      }

      carbo_balance sugar_model_out = sugar_model(year, day, climate.TAir[weather_index],
                                                  climate.PAR[weather_index],
                                                  photosynthesis_per_stem,
                                                  common, parameters,
                                                  D00,
                                                  potential_growth.previous_values.sH,
                                                  resp,
                                                  nitrogen_balance,
                                                  nitrogen_change,
                                                  nitrogen_contrast,
                                                  boolsettings.sperling_model,
                                                  tree_alive,
                                                  boolsettings.storage_grows,
                                                  surplus_c,
                                                  repola_values.needle_mass,
                                                  root_mass,
                                                  mycorrhizal_biomass,

                                                  equilibrium_temperature,
                                                  potential_growth,
                                                  sugar_values_for_next_iteration.sugar,
                                                  sugar_values_for_next_iteration.starch,
                                                  sugar_values_for_next_iteration.previous_values);
      // Saved for the next iteration
      sugar_values_for_next_iteration.previous_values = sugar_model_out.previous_values;
      sugar_values_for_next_iteration.sugar = sugar_model_out.sugar;
      sugar_values_for_next_iteration.starch = sugar_model_out.starch;
      sugar_values_for_next_iteration.storage = sugar_model_out.storage;
      sugar_values_for_next_iteration.nitrogen_capacity = sugar_model_out.nitrogen_capacity;
      parameters.sB0 = sugar_values_for_next_iteration.previous_values.sB0;
      tree_alive = sugar_model_out.previous_values.tree_alive;
      equilibrium_temperature.needles = std::log(sugar_values_for_next_iteration.previous_values.As.needles/sugar_values_for_next_iteration.previous_values.Ad.needles)/(sugar_values_for_next_iteration.starch.B - sugar_values_for_next_iteration.starch.B);
      equilibrium_temperature.phloem = std::log(sugar_values_for_next_iteration.previous_values.As.phloem/sugar_values_for_next_iteration.previous_values.Ad.phloem)/(sugar_values_for_next_iteration.starch.B - sugar_values_for_next_iteration.starch.B);
      equilibrium_temperature.xylem_sh = std::log(sugar_values_for_next_iteration.previous_values.As.xylem_sh/sugar_values_for_next_iteration.previous_values.Ad.xylem_sh)/(sugar_values_for_next_iteration.starch.B - sugar_values_for_next_iteration.starch.B);
      equilibrium_temperature.xylem_st = std::log(sugar_values_for_next_iteration.previous_values.As.xylem_st/sugar_values_for_next_iteration.previous_values.Ad.xylem_st)/(sugar_values_for_next_iteration.starch.B - sugar_values_for_next_iteration.starch.B);
      equilibrium_temperature.roots = std::log(sugar_values_for_next_iteration.previous_values.As.roots/sugar_values_for_next_iteration.previous_values.Ad.roots)/(sugar_values_for_next_iteration.starch.B - sugar_values_for_next_iteration.starch.B);
      if (nitrogen_change) {
        nitrogen_balance = sugar_model_out.nitrogen_balance;

      }

      // TODO: add other organs

      if (final_year%2==0) {
        sugar_values_output.nitrogen_balance.push_back(nitrogen_balance);
        sugar_values_output.nitrogen_capacity_needles.push_back(sugar_model_out.nitrogen_capacity.needles);
        sugar_values_output.nitrogen_capacity_wall.push_back(sugar_model_out.nitrogen_capacity.wall);
        sugar_values_output.nitrogen_capacity_height.push_back(sugar_model_out.nitrogen_capacity.height);
        sugar_values_output.nitrogen_capacity_bud.push_back(sugar_model_out.nitrogen_capacity.bud);
        sugar_values_output.nitrogen_capacity_roots.push_back(sugar_model_out.nitrogen_capacity.roots);

        sugar_values_output.sugar.push_back(sugar_values_for_next_iteration.sugar.needles +
          sugar_values_for_next_iteration.sugar.phloem +
          sugar_values_for_next_iteration.sugar.xylem_sh +
          sugar_values_for_next_iteration.sugar.xylem_st +
          sugar_values_for_next_iteration.sugar.roots);
        sugar_values_output.starch.push_back(sugar_values_for_next_iteration.starch.needles +
          sugar_values_for_next_iteration.starch.phloem +
          sugar_values_for_next_iteration.starch.xylem_sh +
          sugar_values_for_next_iteration.starch.xylem_st +
          sugar_values_for_next_iteration.starch.roots);
        sugar_values_output.storage.push_back(sugar_values_for_next_iteration.storage.needles +
          sugar_values_for_next_iteration.storage.phloem +
          sugar_values_for_next_iteration.storage.xylem_sh +
          sugar_values_for_next_iteration.storage.xylem_st +
          sugar_values_for_next_iteration.storage.roots);

        sugar_values_output.starch_needles.push_back(sugar_values_for_next_iteration.starch.needles);
        sugar_values_output.starch_phloem.push_back(sugar_values_for_next_iteration.starch.phloem);
        sugar_values_output.starch_xylem_sh.push_back(sugar_values_for_next_iteration.starch.xylem_sh);
        sugar_values_output.starch_xylem_st.push_back(sugar_values_for_next_iteration.starch.xylem_st);
        sugar_values_output.starch_roots.push_back(sugar_values_for_next_iteration.starch.roots);

        sugar_values_output.sugar_needles.push_back(sugar_values_for_next_iteration.sugar.needles);
        sugar_values_output.sugar_phloem.push_back(sugar_values_for_next_iteration.sugar.phloem);
        sugar_values_output.sugar_xylem_sh.push_back(sugar_values_for_next_iteration.sugar.xylem_sh);
        sugar_values_output.sugar_xylem_st.push_back(sugar_values_for_next_iteration.sugar.xylem_st);
        sugar_values_output.sugar_roots.push_back(sugar_values_for_next_iteration.sugar.roots);
        sugar_values_output.sugar_mycorrhiza.push_back(sugar_values_for_next_iteration.sugar.mycorrhiza);
        sugar_values_output.sugar_surplus.push_back(sugar_values_for_next_iteration.sugar.surplus);

        sugar_values_output.storage_needles.push_back(sugar_values_for_next_iteration.storage.needles);
        sugar_values_output.storage_phloem.push_back(sugar_values_for_next_iteration.storage.phloem);
        sugar_values_output.storage_xylem_sh.push_back(sugar_values_for_next_iteration.storage.xylem_sh);
        sugar_values_output.storage_xylem_st.push_back(sugar_values_for_next_iteration.storage.xylem_st);
        sugar_values_output.storage_roots.push_back(sugar_values_for_next_iteration.storage.roots);

        if (nitrogen_change) {
          uptake_values_output.ectomycorrhizal_transfer.push_back(sugar_model_out.uptake.ectomycorrhizal_transfer);
          uptake_values_output.root_upatke.push_back(sugar_model_out.uptake.root_upatke);
          uptake_values_output.ectomycorrhizal_upatke.push_back(sugar_model_out.uptake.ectomycorrhizal_upatke);
          uptake_values_output.total_uptake.push_back(sugar_model_out.uptake.total_uptake);
        }

        /*
         * Respiration
         */

        respiration_output.maintenance.push_back(sugar_model_out.resp_main);
        respiration_output.growth.push_back(sugar_model_out.resp_growth);
      }

      /*
       * Actual growth
       */

      double n_rows = ratios.form_factor * parameters.h0 / parameters.cell_l_ew * M_PI * parameters.D0 / parameters.cell_d_ew; // TODO: not used anywhere
      // double GD = // g_sD_T *fD * LD; TODO: where do the parameters come from?
      double ew_cells_vector; // TODO: where does this come from?
      double lw_cells_vector;
      double max_ew_cells = std::max(ew_cells_vector, 0.0); // TODO: make this come from the last iteration?
      double max_lw_cells = std::max(lw_cells_vector, 0.0); // TODO: make this come from the last iteration?

      growth_out actual_growth_out = actual_growth(parameters, common,
                                                   sugar_values_for_next_iteration.storage, potential_growth,
                                                   resp,
                                                   boolsettings.sperling_model,
                                                   sugar_model_out.nitrogen_capacity); // TODO: either than the nitrogen capacity or the nitrogen balance
      // TODO: update the parameters like D0 and h0 that need to be updated
      double growth_and_mortality = actual_growth_out.roots * (0.975 - fS_out); // TODO: more sensible value here!

      ring_width_out ring_width = ring_width_generator(day, previous_ring_width, potential_growth.previous_values, parameters, actual_growth_out.GD);
      previous_ring_width = ring_width;

      double xylem_width;
      double days_in_15_years = 365*15;
      if (weather_index > days_in_15_years) {
        xylem_width = actual_growth_output.ring_width[weather_index] - actual_growth_output.ring_width[weather_index - days_in_15_years];
      } else {
        if (year == start_year & weather_index == 0) {
          xylem_edge_estimate = parameters.xylem_start_estimate;
        } else {
          if (final_year%2!=0) {
            xylem_edge_estimate_internal = parameters.xylem_start_estimate + actual_growth_output.ring_width[weather_index];
            xylem_width = actual_growth_output.ring_width[weather_index] - xylem_edge_estimate_internal;
          } else {
            xylem_edge_estimate = parameters.xylem_start_estimate + actual_growth_output.ring_width[weather_index];
            xylem_width = actual_growth_output.ring_width[weather_index] - xylem_edge_estimate;
          }
        }
      }

      /*
       * Culmative growwth and output
       */

      if (final_year%2!=0) {
        /*
         * Actual growth
         */

        actual_growth_output.height.push_back(actual_growth_out.height);
        actual_growth_output.needles.push_back(actual_growth_out.needles);
        actual_growth_output.roots.push_back(actual_growth_out.roots);
        actual_growth_output.diameter.push_back(actual_growth_out.wall);
        actual_growth_output.bud.push_back(actual_growth_out.bud);
        actual_growth_output.ring_width.push_back(ring_width.tot_mm);

        /*
         * Culmative growth
         */

        if (day == 0) {
          if (year == start_year) {
            culm_growth_internal.height.push_back(height_next_year + actual_growth_out.height);
            culm_growth_internal.diameter.push_back(diameter_next_year + 2*ring_width.tot_mm);
            // cell_wall_density_ew kg C m-3 only early wood uses as early and late wood are currently the same for both parameters
            culm_growth_internal.xylem_st.push_back(M_PI * pow(xylem_width/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth_internal.phloem.push_back(M_PI * pow(1.5/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth_internal.diameter_potential.push_back(diameter_potential_next_year + 2*potential_growth.previous_values.pot_mm);
            culm_growth_internal.needles.push_back(actual_growth_out.needles); // TODO: the life of the needles is actually here, add to the LAI section
            culm_growth_internal.roots.push_back(roots_next_year + growth_and_mortality); // TODO: make this an actual parameter
            culm_growth_internal.mycorrhiza.push_back(mycorrhiza_next_year + growth_and_mortality);
          } else {
            culm_growth_internal.height.push_back(culm_growth_internal.height[weather_index-1] + actual_growth_out.height);
            culm_growth_internal.diameter.push_back(diameter_next_year + 2*ring_width.tot_mm); // Ring width is culmative, but need consistent initial condition
            culm_growth_internal.xylem_st.push_back(M_PI * pow(xylem_width/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth_internal.phloem.push_back(M_PI * pow(1.5/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth_internal.diameter_potential.push_back(diameter_potential_next_year + 2*potential_growth.previous_values.pot_mm);
            culm_growth_internal.needles.push_back(actual_growth_out.needles); // TODO: should there be culmative things here? Fix after the max works
            culm_growth_internal.roots.push_back(culm_growth_internal.roots[weather_index-1] + growth_and_mortality);
            culm_growth_internal.mycorrhiza.push_back(culm_growth_internal.mycorrhiza[weather_index-1] + growth_and_mortality);
          }
        } else {
          culm_growth_internal.height.push_back(culm_growth_internal.height[weather_index-1] + actual_growth_out.height);
          culm_growth_internal.diameter.push_back(diameter_next_year + 2*ring_width.tot_mm);
          culm_growth_internal.xylem_st.push_back(M_PI * pow(xylem_width/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
          culm_growth_internal.phloem.push_back(M_PI * pow(1.5/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
          culm_growth_internal.diameter_potential.push_back(diameter_potential_next_year + 2*potential_growth.previous_values.pot_mm);
          culm_growth_internal.needles.push_back(culm_growth_internal.needles[weather_index-1] + actual_growth_out.needles);
          culm_growth_internal.roots.push_back(culm_growth_internal.roots[weather_index-1] + growth_and_mortality);
          culm_growth_internal.mycorrhiza.push_back(culm_growth_internal.mycorrhiza[weather_index-1] + growth_and_mortality);
        }
      } else {
        if (day == 0) {
          if (year == start_year) {
            culm_growth.height.push_back(height_next_year + actual_growth_out.height);
            culm_growth.diameter.push_back(diameter_next_year + 2*ring_width.tot_mm);
            culm_growth.xylem_st.push_back(M_PI * pow(xylem_width/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth.phloem.push_back(M_PI * pow(1.5/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth.diameter_potential.push_back(diameter_potential_next_year + 2*potential_growth.previous_values.pot_mm);
            culm_growth.needles.push_back(actual_growth_out.needles);
            culm_growth.roots.push_back(roots_next_year + growth_and_mortality);
            culm_growth.mycorrhiza.push_back(mycorrhiza_next_year + growth_and_mortality);
          } else {
            culm_growth.height.push_back(culm_growth.height[weather_index-1] + actual_growth_out.height);
            culm_growth.diameter.push_back(diameter_next_year + 2*ring_width.tot_mm);
            culm_growth.xylem_st.push_back(M_PI * pow(xylem_width/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth.phloem.push_back(M_PI * pow(1.5/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
            culm_growth.diameter_potential.push_back(diameter_potential_next_year + 2*potential_growth.previous_values.pot_mm);
            culm_growth.needles.push_back(actual_growth_out.needles);
            culm_growth.roots.push_back(culm_growth.roots[weather_index-1] + growth_and_mortality);
            culm_growth.mycorrhiza.push_back(culm_growth.mycorrhiza[weather_index-1] + growth_and_mortality);
          }
        } else {
          culm_growth.height.push_back(culm_growth.height[weather_index-1] + actual_growth_out.height);
          culm_growth.diameter.push_back(diameter_next_year + 2*ring_width.tot_mm);
          culm_growth.xylem_st.push_back(M_PI * pow(xylem_width/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
          culm_growth.phloem.push_back(M_PI * pow(1.5/1000, 2) * (height_next_year + actual_growth_out.height) * parameters.cell_wall_density_ew);
          culm_growth.diameter_potential.push_back(diameter_potential_next_year + 2*potential_growth.previous_values.pot_mm);
          culm_growth.needles.push_back(culm_growth.needles[weather_index-1] + actual_growth_out.needles);
          culm_growth.roots.push_back(culm_growth.roots[weather_index-1] + growth_and_mortality);
          culm_growth.mycorrhiza.push_back(culm_growth.mycorrhiza[weather_index-1] + growth_and_mortality);
        }

        culm_growth.tree_alive.push_back(tree_alive);
      }

      /*
       * Updating some parameters
       */

      // TODO: what does this do?
      HH = potential_growth.previous_values.HH;
      potenital_growth_use.push_back(potential_growth.use);

      if (final_year%2==0) {
        years.push_back(year);
        days.push_back(day+1);

        if (day == days_per_year-1) {
          diameter_next_year = culm_growth.diameter[weather_index];
          diameter_potential_next_year = culm_growth.diameter_potential[weather_index];
        }
      }

      GPP_sum_yesterday = GPP_sum;
    } // End of the days loop

    // Updating parameters!

    last_year_HH = HH;

    if (final_year%2==0) {
      GPP_previous_sum.push_back(std::accumulate(photosynthesis_output.GPP.begin() + days_gone + 182, photosynthesis_output.GPP.begin() + days_gone + 245 + 1, 0.0));

      // TOOD: should this be every year or just every other?
      max_needles = *std::max_element(culm_growth.needles.begin() + days_gone,
                                      culm_growth.needles.begin() + days_gone + days_per_year);

      days_gone = days_gone + days_per_year;

      last_cohorts.year_1 = needles_cohorts.year_1; // TODO: currently the growth doesn't really have an effect on this - should it?
      last_cohorts.year_2 = needles_cohorts.year_2;
      last_cohorts.year_3 = needles_cohorts.year_3;
      // TODO: make it possible to make this for 5 years?

      parameters.sugar0 = sugar_values_for_next_iteration.sugar.needles + sugar_values_for_next_iteration.sugar.phloem + sugar_values_for_next_iteration.sugar.roots + sugar_values_for_next_iteration.sugar.xylem_sh + sugar_values_for_next_iteration.sugar.xylem_st;
      parameters.sugar_needles0 = sugar_values_for_next_iteration.sugar.needles;
      parameters.sugar_phloem0 = sugar_values_for_next_iteration.sugar.phloem;
      parameters.sugar_roots0 = sugar_values_for_next_iteration.sugar.roots;
      parameters.sugar_xylem_sh0 = sugar_values_for_next_iteration.sugar.xylem_sh;
      parameters.sugar_xylem_st0 = sugar_values_for_next_iteration.sugar.xylem_st;

      parameters.starch0 = sugar_values_for_next_iteration.starch.needles + sugar_values_for_next_iteration.starch.phloem + sugar_values_for_next_iteration.starch.roots + sugar_values_for_next_iteration.starch.xylem_sh + sugar_values_for_next_iteration.starch.xylem_st;
      parameters.starch_needles0 = sugar_values_for_next_iteration.starch.needles;
      parameters.starch_phloem0 = sugar_values_for_next_iteration.starch.phloem;
      parameters.starch_roots0 = sugar_values_for_next_iteration.starch.roots;
      parameters.starch_xylem_sh0 = sugar_values_for_next_iteration.starch.xylem_sh;
      parameters.starch_xylem_st0 = sugar_values_for_next_iteration.starch.xylem_st;
    }

    final_year = final_year + 1;
  } // End of the years loop

  ///////////////////////
  // Output
  ///////////////////////

  // Dataframes

  Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::_["year"] = years,
                                               Rcpp::_["day"] = days,
                                               Rcpp::_["g"] = potential_growth_output.g,
                                               Rcpp::_["bud_growth_potential"] = potential_growth_output.bud,
                                               Rcpp::_["wall_growth_potential"] = potential_growth_output.diameter,
                                               Rcpp::_["needle_growth_potential"] = potential_growth_output.needles,
                                               Rcpp::_["roots_growth_potential"] = potential_growth_output.roots,
                                               Rcpp::_["height_growth_potential"] = potential_growth_output.height,
                                               Rcpp::_["bud_growth"] = actual_growth_output.bud,
                                               Rcpp::_["diameter_growth"] = actual_growth_output.diameter,
                                               Rcpp::_["needle_growth"] = actual_growth_output.needles,
                                               Rcpp::_["root_growth"] = actual_growth_output.roots,
                                               Rcpp::_["height_growth"] = actual_growth_output.height,
                                               Rcpp::_["respiration_growth"] = respiration_output.growth,
                                               Rcpp::_["respiration_maintenance"] = respiration_output.maintenance,
                                               Rcpp::_["ring_width"] = actual_growth_output.ring_width,
                                               Rcpp::_["mm_potenital"] = potential_growth_output.ring_width);

  Rcpp::DataFrame df2 = Rcpp::DataFrame::create(Rcpp::_["sugar"] = sugar_values_output.sugar,
                                                Rcpp::_["starch"] = sugar_values_output.starch,
                                                Rcpp::_["storage"] = sugar_values_output.storage,
                                                Rcpp::_["starch_needles"] = sugar_values_output.starch_needles,
                                                Rcpp::_["starch_phloem"] = sugar_values_output.starch_phloem,
                                                Rcpp::_["starch_xylem_sh"] = sugar_values_output.starch_xylem_sh,
                                                Rcpp::_["starch_xylem_st"] = sugar_values_output.starch_xylem_st,
                                                Rcpp::_["starch_roots"] = sugar_values_output.starch_roots,
                                                Rcpp::_["sugar_needles"] = sugar_values_output.sugar_needles,
                                                Rcpp::_["sugar_phloem"] = sugar_values_output.sugar_phloem,
                                                Rcpp::_["sugar_xylem_sh"] = sugar_values_output.sugar_xylem_sh,
                                                Rcpp::_["sugar_xylem_st"] = sugar_values_output.sugar_xylem_st,
                                                Rcpp::_["sugar_roots"] = sugar_values_output.sugar_roots,
                                                Rcpp::_["sugar_to_mycorrhiza"] = sugar_values_output.sugar_mycorrhiza,
                                                Rcpp::_["sugar_out_surplus"] = sugar_values_output.sugar_surplus,
                                                Rcpp::_["storage_term_needles"] = sugar_values_output.storage_needles,
                                                Rcpp::_["storage_term_phloem"] = sugar_values_output.storage_phloem,
                                                Rcpp::_["storage_term_xylem_sh"] = sugar_values_output.storage_xylem_sh,
                                                Rcpp::_["storage_term_xylem_st"] = sugar_values_output.storage_xylem_st,
                                                Rcpp::_["storage_term_roots"] = sugar_values_output.storage_roots,
                                                Rcpp::_["n_E_pot"] = growth_values_for_next_iteration.n_E_pot,
                                                Rcpp::_["n_W_pot"] = growth_values_for_next_iteration.n_W_pot,
                                                Rcpp::_["n_M_pot"] = growth_values_for_next_iteration.n_M_pot);

  Rcpp::DataFrame df3 = Rcpp::DataFrame::create(Rcpp::_["GPP"] = photosynthesis_output.GPP,
                                                Rcpp::_["ET"] = photosynthesis_output.ET,
                                                Rcpp::_["SoilWater"] = photosynthesis_output.SoilWater,
                                                Rcpp::_["fAPAR"] = photosynthesis_output.fAPAR);

  Rcpp::DataFrame df4 = Rcpp::DataFrame::create(Rcpp::_["culm_growth_height"] = culm_growth.height,
                                                Rcpp::_["culm_growth_needles"] = culm_growth.needles,
                                                Rcpp::_["culm_growth_diameter"] = culm_growth.diameter,
                                                Rcpp::_["culm_growth_diameter_potential"] = culm_growth.diameter_potential,
                                                Rcpp::_["culm_growth_roots"] = culm_growth.roots,
                                                Rcpp::_["tree_alive"] = culm_growth.tree_alive,
                                                Rcpp::_["LAI"] = LAI_within_year_vector,
                                                Rcpp::_["nitrogen_balance"] = sugar_values_output.nitrogen_balance,
                                                Rcpp::_["nitrogen_capacity_needles"] = sugar_values_output.nitrogen_capacity_needles,
                                                Rcpp::_["nitrogen_capacity_wall"] = sugar_values_output.nitrogen_capacity_wall,
                                                Rcpp::_["nitrogen_capacity_height"] = sugar_values_output.nitrogen_capacity_height,
                                                Rcpp::_["nitrogen_capacity_bud"] = sugar_values_output.nitrogen_capacity_bud,
                                                Rcpp::_["nitrogen_capacity_roots"] = sugar_values_output.nitrogen_capacity_roots);

  if (nitrogen_change) {
    Rcpp::DataFrame df5 = Rcpp::DataFrame::create(Rcpp::_["root_upatke"] = uptake_values_output.root_upatke,
                                                  Rcpp::_["ectomycorrhizal_upatke"] = uptake_values_output.ectomycorrhizal_upatke,
                                                  Rcpp::_["total_uptake"] = uptake_values_output.total_uptake,
                                                  Rcpp::_["ectomycorrhizal_transfer"] = uptake_values_output.ectomycorrhizal_transfer);

    return Rcpp::List::create(Rcpp::_["Growth"] = df,
                              Rcpp::_["Sugar"] = df2,
                              Rcpp::_["Preles"] = df3,
                              Rcpp::_["Culm_Growth"] = df4,
                              Rcpp::_["Uptake"] = df5);
  } else {

    return Rcpp::List::create(Rcpp::_["Growth"] = df,
                              Rcpp::_["Sugar"] = df2,
                              Rcpp::_["Preles"] = df3,
                              Rcpp::_["Culm_Growth"] = df4);
  }

}
