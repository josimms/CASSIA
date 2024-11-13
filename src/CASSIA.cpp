#include "CASSIA.h"

int leap_year(int year)
{
  if (year % 4 != 0) {
    return 365;
  }
  else {
    return 366;
  }
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

  CASSIA_parameter_test(parameters);
  // CASSIA_common_test(pCASSIA_common);
  // CASSIA_ratios_test(pCASSIA_ratios, "Hyde");

  Settings boolsettings = parseSettings(settings);

  /*
   * Weather input made into vectors
   */

  std::vector<double> PAR = weather["PAR"];
  std::vector<double> VPD = weather["VPD"];
  std::vector<double> fAPAR = weather["fAPAR"];
  std::vector<double> Nitrogen = weather["Nitrogen"];
  std::vector<double> Photosynthesis_IN = weather["P"];
  std::vector<double> TAir = weather["T"];
  std::vector<double> TSoil_A = weather["TSA"];
  std::vector<double> TSoil_B = weather["TSB"];
  std::vector<double> Soil_Moisture = weather["MB"];
  std::vector<double> Precip  = weather["Rain"];
  std::vector<double> CO2 = weather["CO2"];

  /*
   * Structures set up
   */

  // Forward init values (previous day values) as first values of result vectors

  bool tree_alive = TRUE;

  // std::vector<double> SW, Canopywater, SOG, S;

  /*
   * SW.push_back(SWinit);
   Canopywater.push_back(CWinit);
   SOG.push_back(SOGinit);
   S.push_back(Sinit);
   */

  double CH = parameters.density_tree * parameters.carbon_share;
  double M_suc = 12 * common.M_C + 22 * common.M_H + 11 * common.M_O;

  carbo_tracker Ad;
  Ad.needles = parameters.Ad0_needles;
  Ad.phloem = parameters.Ad0_phloem;
  Ad.xylem_sh = parameters.Ad0_xylem_sh;
  Ad.xylem_st = parameters.Ad0_xylem_st;
  Ad.roots = parameters.Ad0_roots;

  repola_out repola_values;
  if (needle_mass_in == 0) { // The value of this should be 0 if you want the needle value to be calculated
    repola_values = repola(parameters); // TODO: fix
  } else {
    repola_values.needle_mass = needle_mass_in;
  }

  /*
   * Vectors between iterations
   */

  growth_values_out growth_values_for_next_iteration;
  carbo_balance sugar_values_for_next_iteration;
  ring_width_out previous_ring_width;
  carbo_balance original_parameters;

  /*
   * Vectors for the outputs
   */

  growth_vector potential_growth_output;
  growth_vector actual_growth_output;
  sugar_values_vector sugar_values_output;
  photo_out_vector photosynthesis_output;
  resp_vector respiration_output;
  needle_cohorts last_cohorts;
  double last_year_HH;
  double last_year_maxN;
  double GPP_mean;
  std::vector<double> GPP_previous_sum;
  GPP_previous_sum.push_back(481.3); // TODO; make this a variable input, rather than this 2015 value
  double respiration_maintanence;
  std::vector<double> potenital_growth_use;

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
  for (int year : years_for_runs)  {

    // std::cout << " Year " << year;

    /*
     * Daily output
     */

    carbo_tracker carbo_tracker_vector;
    xylogensis_out xylogensis_vector;
    double fS;

    /*
     * Yearly initialization
     */

    // Temperature equilibrium for the sugar model
    //	# Compute initial Te by the mean temperature for the first week of # Semptemver plus 3C (for the exponential nature of the curves), original was October
    double equilibrium_temperature = (TAir[244] + TAir[245] + TAir[246] + TAir[247] + TAir[248] + TAir[249] + TAir[250] + TAir[251]) / 7 + 3;

    // B0, D00 and h00
    double B0 = M_PI/4 * pow(parameters.D0, 2);
    double D00 = parameters.D0;
    double h00 = parameters.h0;
    if (boolsettings.xylogensis_option) {
      double LH0 = parameters.h_increment / (0.5 * parameters.sHc);
      double LN0 = parameters.n_length / (0.5 * parameters.sNc);
      double LR0 = 2 * parameters.m_R_tot / parameters.sRc; // TODO: check if these parameters are updated somewhere
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
      needles_cohorts.year_1 = 31.24535 / parameters.n_length * repola_values.needle_mass / 3;
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
    if (boolsettings.CASSIA_graphs) {
      days_per_year = 365;
    }
    /*
     * Yearly initial conditions updated
     */

    // Set up the initial conditions
    yearly_in yearly = yearly_initial_conditions(days_per_year); // TODO: growth should be added to this!

    /*
     * DAYS LOOP
     */
    for (int day = 0; day < days_per_year; day++) {

      /*
       * PHOTOSYNTHESIS
       *
       * There are yearly, but not daily dependencies other than environmental states here!
       */

      // Uses the method from Tian 2021
      // Extinction coefficient 0.52 is from Tian 2021 as well
      // LAI value is fairly constant if we look at Rautiainen 2012, LAI ~ 3
      double LAI = 3; // TODO
      double f_modifer = needles_last/growth_values_for_next_iteration.max_N;
      double LAI_within_year;
      if (day < 182) { // TODO: I decided that the start of July is the end of spring
        LAI_within_year = 2.0/3.0*LAI + f_modifer*1.0/3.0*LAI;
      } else if (day > 244) { // TODO: I decided that the end of august is the start of autumn
        LAI_within_year = 2.0/3.0*LAI + fS*1.0/3.0*LAI;
      } else {
        LAI_within_year = LAI;
      }
      double fAPAR = (1 - std::exp(-0.52 * LAI_within_year));  // TODO: Check this is sensible

      double photosynthesis_per_stem;
      photosynthesis_out photosynthesis;
      if (boolsettings.photosynthesis_as_input) {
        photosynthesis.GPP = Photosynthesis_IN[day];
        photosynthesis_per_stem = Photosynthesis_IN[day] / 1010 * 10000/1000;
      } else {
        photosynthesis.GPP = 0;
        photosynthesis.ET = 0;
        photosynthesis.fS = 0;
        photosynthesis.SoilWater = 0;
        // call_preles(days_per_year, day,
        //           PAR[count], TAir[count], VPD[count], Precip[count],
        //           CO2[count], fAPAR, Nitrogen[count],
        //           parSite, parGPP, parET, parSnowRain,
        //           parWater, parN, etmodel);
        // TODO: 5.6...
        photosynthesis.GPP = 5.6 * photosynthesis.GPP; // g C m-2 per day, so no conversion is needed!
        photosynthesis_per_stem = photosynthesis.GPP / 1010 * 10000/1000;
        fS = photosynthesis.fS;
      }
      if (day == 0) {
        GPP_sum = 0.0;
      } else if (day <= 182) {
        GPP_sum = GPP_sum_yesterday;
      } else if (day > 182 && day <= 244) {
        GPP_sum = GPP_sum_yesterday + photosynthesis.GPP;
      } else if (day > 245) {
        GPP_sum = GPP_sum_yesterday;
      }

      std::cout << "photosynthesis.GPP: " << photosynthesis.GPP << " Photosynthesis: " << Photosynthesis_IN[day] << "\n";

      /*
       * Potential Growth
       *
       * In terms of the adaptation from the R code, the potential values are not altered by daily processes so still calculate them for a year
       */

      double en_pot_growth_old;
      if (day < 11) {
        en_pot_growth_old = 0.0;
      } else {
        en_pot_growth_old = potenital_growth_use[day-11];
      }

      GPP_mean = 463.8833; // TODO: move when I understand GPP_sum
      growth_out potential_growth = growth(day, year, TAir[day], TSoil_A[day], TSoil_B[day], Soil_Moisture[day], photosynthesis.GPP, GPP_ref[day],
                                           boolsettings.root_as_Ding, boolsettings.xylogensis_option, boolsettings.environmental_effect_xylogenesis, boolsettings.sD_estim_T_count,
                                           common, parameters, ratios,
                                           CH, B0, en_pot_growth_old, GPP_mean, GPP_previous_sum[year-start_year],
                                                                                                boolsettings.LH_estim, boolsettings.LN_estim, boolsettings.LD_estim,
                                                                                                // Last iteration value
                                                                                                growth_values_for_next_iteration, last_year_HH,
                                                                                                days_per_year);
      // Saved for the next iteration
      growth_values_for_next_iteration = potential_growth.previous_values;

      /*
       * Respiration
       *
       * There are yearly, but not daily dependencies other than weather conditions here!
       */

      respiration_out resp = respiration(day, parameters, ratios, repola_values,
                                         TAir[day], TSoil_A[day],
                                                             boolsettings.temp_rise, boolsettings.Rm_acclimation, boolsettings.mN_varies,
                                                             // parameters that I am not sure about
                                                             B0);

      /*
       * Sugar
       */

      // TODO; need to check the indexes!
      if ((day == 0) & (year == start_year)) {
        sugar_values_for_next_iteration.sugar.needles = original_parameters.sugar.needles = parameters.sugar_needles0;
        sugar_values_for_next_iteration.sugar.phloem = original_parameters.sugar.phloem = parameters.sugar_phloem0;
        sugar_values_for_next_iteration.sugar.roots = original_parameters.sugar.roots = parameters.sugar_roots0;
        sugar_values_for_next_iteration.sugar.xylem_sh = original_parameters.sugar.xylem_sh = parameters.sugar_xylem_sh0;
        sugar_values_for_next_iteration.sugar.xylem_st = original_parameters.sugar.xylem_st = parameters.sugar_xylem_st0;
        sugar_values_for_next_iteration.sugar.mycorrhiza = original_parameters.sugar.mycorrhiza = 0; // TODO: think about this

        sugar_values_for_next_iteration.starch.needles = original_parameters.starch.needles = parameters.starch_needles0;
        sugar_values_for_next_iteration.starch.phloem = original_parameters.starch.phloem = parameters.starch_phloem0;
        sugar_values_for_next_iteration.starch.roots = original_parameters.starch.roots = parameters.starch_roots0;
        sugar_values_for_next_iteration.starch.xylem_sh = original_parameters.starch.xylem_sh = parameters.starch_xylem_sh0;
        sugar_values_for_next_iteration.starch.xylem_st = original_parameters.starch.xylem_st = parameters.starch_xylem_st0;
        sugar_values_for_next_iteration.starch.mycorrhiza = original_parameters.starch.mycorrhiza = 0;
      } else if ((day == 0) & (year != start_year)) {
        sugar_values_for_next_iteration.sugar.needles = parameters.sugar_needles0;
        sugar_values_for_next_iteration.sugar.phloem = parameters.sugar_phloem0;
        sugar_values_for_next_iteration.sugar.roots = parameters.sugar_roots0;
        sugar_values_for_next_iteration.sugar.xylem_sh = parameters.sugar_xylem_sh0;
        sugar_values_for_next_iteration.sugar.xylem_st = parameters.sugar_xylem_st0;
        sugar_values_for_next_iteration.sugar.mycorrhiza = 0;

        sugar_values_for_next_iteration.starch.needles = parameters.starch_needles0;
        sugar_values_for_next_iteration.starch.phloem = parameters.starch_phloem0;
        sugar_values_for_next_iteration.starch.roots = parameters.starch_roots0;
        sugar_values_for_next_iteration.starch.xylem_sh = parameters.starch_xylem_sh0;
        sugar_values_for_next_iteration.starch.xylem_st = parameters.starch_xylem_st0;
        sugar_values_for_next_iteration.starch.mycorrhiza = 0;
      } else {
        carbo_balance sugar_model_out = sugar_model(day, TAir[day], photosynthesis_per_stem,
                                                    common, parameters,
                                                    D00,
                                                    potential_growth.previous_values.sH,
                                                    resp,
                                                    boolsettings.sperling_model,
                                                    tree_alive,
                                                    boolsettings.storage_grows,
                                                    repola_values.needle_mass,
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
        parameters.sB0 = sugar_values_for_next_iteration.previous_values.sB0;
        tree_alive = sugar_model_out.previous_values.tree_alive;
      }

      /*
       if (sugar_values_for_next_iteration.sugar.needles < 0 || sugar_values_for_next_iteration.sugar.phloem < 0 ||
       sugar_values_for_next_iteration.sugar.xylem_sh < 0 || sugar_values_for_next_iteration.sugar.xylem_st < 0 ||
       sugar_values_for_next_iteration.sugar.roots < 0) {
       std::cout << "SUGAR IS NEGATIVE, STOP SIMULATION!\n";
       return(Rcpp::DataFrame::create(0));
       }
       */

      /*
       * Actual growth
       */

      double n_rows = ratios.form_factor * parameters.h0 / parameters.cell_l_ew * M_PI * parameters.D0 / parameters.cell_d_ew;
      // double GD = // g_sD_T *fD * LD; TODO: where do the parameters come from?
      double ew_cells_vector; // TODO: where does this come from?
      double lw_cells_vector;
      double max_ew_cells = std::max(ew_cells_vector, 0.0); // TODO: make this come from the last iteration?
      double max_lw_cells = std::max(lw_cells_vector, 0.0); // TODO: make this come from the last iteration?

      // TODO: storage
      growth_out actual_growth_out = actual_growth(parameters, common,
                                                   sugar_values_for_next_iteration.storage, potential_growth,
                                                   resp,
                                                   boolsettings.sperling_model);
      // TODO: update the parameters like D0 and h0 that need to be updated

      ring_width_out ring_width = ring_width_generator(day, previous_ring_width, potential_growth.previous_values, parameters, actual_growth_out.GD);
      previous_ring_width = ring_width;

      /*
       * Output
       */
      HH = potential_growth.previous_values.HH;
      needles_last = potential_growth.needles;
      potenital_growth_use.push_back(potential_growth.use);
      if (!tree_alive) {
        // std::cout << "The tree is dead due to sugar storage\n";
      }

      GPP_sum_yesterday = GPP_sum;

      if (final_year%2==0) {
        years.push_back(year);
        days.push_back(day+1);

        respiration_output.maintenance.push_back(actual_growth_out.respiration_maintenance);
        respiration_output.growth.push_back(actual_growth_out.respiration_growth);

        potential_growth_output.height.push_back(potential_growth.height);
        potential_growth_output.needles.push_back(potential_growth.needles);
        potential_growth_output.roots.push_back(potential_growth.roots);
        potential_growth_output.diameter.push_back(potential_growth.diameter);
        potential_growth_output.bud.push_back(potential_growth.bud);
        potential_growth_output.g.push_back(potential_growth.g);
        potential_growth_output.en_pot_growth.push_back(potential_growth.use);

        actual_growth_output.height.push_back(actual_growth_out.height);
        actual_growth_output.needles.push_back(actual_growth_out.needles);
        actual_growth_output.roots.push_back(actual_growth_out.roots);
        actual_growth_output.diameter.push_back(actual_growth_out.wall);
        actual_growth_output.bud.push_back(actual_growth_out.bud);
        actual_growth_output.ring_width.push_back(ring_width.tot_mm);

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
        sugar_values_output.starch_needles.push_back(sugar_values_for_next_iteration.starch.needles);
        sugar_values_output.starch_phloem.push_back(sugar_values_for_next_iteration.starch.phloem);
        sugar_values_output.starch_xylem_sh.push_back(sugar_values_for_next_iteration.starch.xylem_sh);
        sugar_values_output.starch_xylem_st.push_back(sugar_values_for_next_iteration.starch.xylem_st);
        sugar_values_output.starch_roots.push_back(sugar_values_for_next_iteration.starch.roots);
        sugar_values_output.starch_mycorrhiza.push_back(sugar_values_for_next_iteration.starch.mycorrhiza);
        sugar_values_output.sugar_needles.push_back(sugar_values_for_next_iteration.sugar.needles);
        sugar_values_output.sugar_phloem.push_back(sugar_values_for_next_iteration.sugar.phloem);
        sugar_values_output.sugar_xylem_sh.push_back(sugar_values_for_next_iteration.sugar.xylem_sh);
        sugar_values_output.sugar_xylem_st.push_back(sugar_values_for_next_iteration.sugar.xylem_st);
        sugar_values_output.sugar_roots.push_back(sugar_values_for_next_iteration.sugar.roots);
        sugar_values_output.sugar_mycorrhiza.push_back(sugar_values_for_next_iteration.sugar.mycorrhiza);

        photosynthesis_output.GPP.push_back(photosynthesis.GPP);
        photosynthesis_output.ET.push_back(photosynthesis.ET);
        photosynthesis_output.SoilWater.push_back(photosynthesis.SoilWater);
        photosynthesis_output.fS.push_back(photosynthesis.fS);
      }
    }

    // Updating parameters!

    last_year_HH = HH;

    if (final_year%2==0) {
      GPP_mean = 463.8833; // TODO: should change this!
      GPP_previous_sum.push_back(GPP_sum);

      last_cohorts.year_1 = needles_cohorts.year_1; // TODO: currently the growth doesn't really have an effect on this - should it?
      last_cohorts.year_2 = needles_cohorts.year_2;
      last_cohorts.year_3 = needles_cohorts.year_3;
      // TODO: make it possible to make this for 5 years?

      parameters.sugar_needles0 = sugar_values_for_next_iteration.sugar.needles;
      parameters.sugar_phloem0 = sugar_values_for_next_iteration.sugar.phloem;
      parameters.sugar_roots0 = sugar_values_for_next_iteration.sugar.roots;
      parameters.sugar_xylem_sh0 = sugar_values_for_next_iteration.sugar.xylem_sh;
      parameters.sugar_xylem_st0 = sugar_values_for_next_iteration.sugar.xylem_st;

      parameters.starch_needles0 = sugar_values_for_next_iteration.starch.needles;
      parameters.starch_phloem0 = sugar_values_for_next_iteration.starch.phloem;
      parameters.starch_roots0 = sugar_values_for_next_iteration.starch.roots;
      parameters.starch_xylem_sh0 = sugar_values_for_next_iteration.starch.xylem_sh;
      parameters.starch_xylem_st0 = sugar_values_for_next_iteration.starch.xylem_st;
    }

    // TODO: need to add the growth of things here!
    final_year = final_year + 1;

    if (year == final_year + 1) {
      parameters.sugar_needles0 = original_parameters.sugar.needles;
      parameters.sugar_phloem0 = original_parameters.sugar.phloem;
      parameters.sugar_roots0 = original_parameters.sugar.roots;
      parameters.sugar_xylem_sh0 = original_parameters.sugar.xylem_sh;
      parameters.sugar_xylem_st0 = original_parameters.sugar.xylem_st;

      parameters.starch_needles0 = original_parameters.starch.needles;
      parameters.starch_phloem0 = original_parameters.starch.phloem;
      parameters.starch_roots0 = original_parameters.starch.roots;
      parameters.starch_xylem_sh0 = original_parameters.starch.xylem_sh;
      parameters.starch_xylem_st0 = original_parameters.starch.xylem_st;
    }
  }

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
                                               Rcpp::_["ring_width"] = actual_growth_output.ring_width);
  Rcpp::DataFrame df2 = Rcpp::DataFrame::create(Rcpp::_["sugar"] = sugar_values_output.sugar,
                                                Rcpp::_["starch"] = sugar_values_output.starch,
                                                Rcpp::_["starch_needles"] = sugar_values_output.starch_needles,
                                                Rcpp::_["starch_phloem"] = sugar_values_output.starch_phloem,
                                                Rcpp::_["starch_xylem_sh"] = sugar_values_output.starch_xylem_sh,
                                                Rcpp::_["starch_xylem_st"] = sugar_values_output.starch_xylem_st,
                                                Rcpp::_["starch_roots"] = sugar_values_output.starch_roots,
                                                Rcpp::_["starch_mycorrhiza"] = sugar_values_output.starch_mycorrhiza,
                                                Rcpp::_["sugar_needles"] = sugar_values_output.sugar_needles,
                                                Rcpp::_["sugar_phloem"] = sugar_values_output.sugar_phloem,
                                                Rcpp::_["sugar_xylem_sh"] = sugar_values_output.sugar_xylem_sh,
                                                Rcpp::_["sugar_xylem_st"] = sugar_values_output.sugar_xylem_st,
                                                Rcpp::_["sugar_roots"] = sugar_values_output.sugar_roots,
                                                Rcpp::_["sugar_mycorrhiza"] = sugar_values_output.sugar_mycorrhiza,
                                                Rcpp::_["n_E_pot"] = growth_values_for_next_iteration.n_E_pot,
                                                Rcpp::_["n_W_pot"] = growth_values_for_next_iteration.n_W_pot,
                                                Rcpp::_["n_M_pot"] = growth_values_for_next_iteration.n_M_pot);
  Rcpp::DataFrame df3 = Rcpp::DataFrame::create(Rcpp::_["GPP"] = photosynthesis_output.GPP,
                                                Rcpp::_["ET"] = photosynthesis_output.ET,
                                                Rcpp::_["SoilWater"] = photosynthesis_output.SoilWater);

  return Rcpp::List::create(Rcpp::_["Growth"] = df,
                            Rcpp::_["Sugar"] = df2,
                            Rcpp::_["Preles"] = df3);

}
