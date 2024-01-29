#include "CASSIA.h"

carbo_tracker storage_carbohydrate(carbo_tracker critical_W, CASSIA_parameters parameters) {
  carbo_tracker out;

  out.needles = 1 / (1 - 1/std::exp(parameters.alfa_needles) * (critical_W.needles - parameters.Wala_needles));
  out.phloem = 1 / (1 - 1/std::exp(parameters.alfa_phloem) * (critical_W.phloem - parameters.Wala_phloem));
  out.xylem_sh = 1 / (1 - 1/std::exp(parameters.alfa_xylem_sh) * (critical_W.xylem_sh - parameters.Wala_xylem_sh));
  out.xylem_st = 1 / (1 - 1/std::exp(parameters.alfa_xylem_st) * (critical_W.xylem_st - parameters.Wala_xylem_st));
  out.roots = 1 / (1 - 1/std::exp(parameters.alfa_roots) * (critical_W.roots - parameters.Wala_roots));

  return out;
}

double storage_update(double alfa, double sugar, double starch, double Wala, bool tree_alive) {
  double out;
  // The checks here make sense, as the model should compensate for the sugar lost at the end of each iteration rather than the beginning
  // Therefore, there should always be a positive value of sugar at the beginning of an iteration (and hopefully generally)
  if (std::isnan(sugar)) {
    // std::cout << "Sugar is NaN ";
    out = 0;
  } else if (std::isnan(starch)) {
    // std::cout << "Starch is NaN ";
    out = 0;
  } else if (sugar < 0) {
    // std::cout << "Sugar is negative ";
    out = 0;
  } else if (starch < 0) {
    // std::cout << "Starch is negative ";
    out = 0;
  } else if (!tree_alive) {
    out = 0;
  } else {
    double ak = 1 / (1 - 1/std::exp(alfa * (sugar + starch - Wala)));
    double comparison = std::min(1.0, ak * (1 - 1/std::exp(alfa * (sugar + starch - Wala))));
    out = std::max(0.0, comparison);
  }
  return out;
}

carbo_tracker As_initiliser(carbo_tracker Ad, double equilibrium_temperature, double Bd, double Bs)
{
  carbo_tracker As;
  As.needles = Ad.needles * std::exp(equilibrium_temperature * (Bd - Bs));
  As.phloem = Ad.phloem * std::exp(equilibrium_temperature * (Bd - Bs));
  As.xylem_sh = Ad.xylem_sh * std::exp(equilibrium_temperature * (Bd - Bs));
  As.xylem_st = Ad.xylem_st * std::exp(equilibrium_temperature * (Bd - Bs));
  As.roots = Ad.roots * std::exp(equilibrium_temperature * (Bd - Bs));

  return As;
}

double emergancy(double sugar, double starch, double tau_emergancy, double lower_bound) {
  double out;
  if (starch > 0) {
    if (sugar < lower_bound) {
      double comaprison = std::max((lower_bound - sugar) / tau_emergancy, 0.0);
      out = std::min(starch, comaprison);
    } else {
      out = 0;
    }
  } else {
    out = 0;
  }

  return out;
}

/*
 * Some of the parameters are site dependent, could define the site side in the main function code so I just import the right variables into this code
 *
 * needles mass per year, should define the year before this function!
 */

carbo_balance sugar_model(int day,
                          double TAir,
                          double PF,

                          CASSIA_common common,
                          CASSIA_parameters parameters,

                          double D00,
                          double sH,
                          respiration_out resp,

                          bool sperling_sugar_model,
                          bool tree_alive,
                          bool storage_grows,
                          double needles_mass, // Repola
                          double temperature_equilibrium,

                          growth_out pot_growth,

                          carbo_tracker sugar,
                          carbo_tracker starch,

                          carbo_values_out parameters_in) {


  carbo_tracker storage_term;
  // TODO: consider these terms in the calibration
  storage_term.needles = storage_update(parameters.alfa_needles, sugar.needles, starch.needles, parameters.lower_bound_needles, tree_alive);
  storage_term.phloem = storage_update(parameters.alfa_phloem, sugar.phloem, starch.phloem, parameters.lower_bound_phloem, tree_alive);
  storage_term.roots = storage_update(parameters.alfa_roots, sugar.roots, starch.roots, parameters.lower_bound_roots, tree_alive);
  storage_term.xylem_sh = storage_update(parameters.alfa_xylem_sh, sugar.xylem_sh, starch.xylem_sh, parameters.lower_bound_xylem_sh, tree_alive);
  storage_term.xylem_st = storage_update(parameters.alfa_xylem_st, sugar.xylem_st, starch.xylem_st, parameters.lower_bound_xylem_st, tree_alive);

  double sB0;
  carbo_tracker As_new, Ad_new;

  if (sperling_sugar_model) {
    if (storage_grows) {
      parameters.lower_bound_needles = parameters.HN0 / D00 * parameters.lower_bound_needles;
      parameters.lower_bound_phloem =  parameters.D0 / D00 * parameters.lower_bound_phloem;
      parameters.lower_bound_roots =  parameters.LR0 / D00 * parameters.lower_bound_roots;
      parameters.lower_bound_xylem_sh =  parameters.LR0 / D00 * parameters.lower_bound_xylem_sh;
      parameters.lower_bound_xylem_st =  parameters.LR0 / D00 * parameters.lower_bound_xylem_st;

      // critical_W.needles = parameters.HN0 / D00 * critical_W.needles;
      // critical_W.phloem =  parameters.D0 / D00 * critical_W.phloem;
      // critical_W.roots =  parameters.LR0 / D00 * critical_W.roots;
      // critical_W.xylem_sh =  parameters.LR0 / D00 * critical_W.xylem_sh;
      // critical_W.xylem_st =  parameters.LR0 / D00 * critical_W.xylem_st;
    }

    /*
     * PARAMETER UPDATES!
     */

    carbo_tracker Kd;
    carbo_tracker Ks;
    // std::cout << " starch.B " << starch.B;
    Kd.needles = parameters_in.Ad.needles*std::exp(starch.B*TAir);
    Ks.needles = parameters_in.As.needles*std::exp(sugar.B*TAir);
    Kd.phloem = parameters_in.Ad.phloem*std::exp(starch.B*TAir);
    Ks.phloem = parameters_in.As.phloem*std::exp(sugar.B*TAir);
    Kd.roots = parameters_in.Ad.roots*std::exp(starch.B*TAir);
    Ks.roots = parameters_in.As.roots*std::exp(sugar.B*TAir);
    Kd.xylem_sh = parameters_in.Ad.xylem_sh*std::exp(starch.B*TAir);
    Ks.xylem_sh = parameters_in.As.xylem_sh*std::exp(sugar.B*TAir);
    Kd.xylem_st = parameters_in.Ad.xylem_st*std::exp(starch.B*TAir);
    Ks.xylem_st = parameters_in.As.xylem_st*std::exp(sugar.B*TAir);

    /*
     * The differences are normalised by a multiplier which represents the average difference in magnitude between the two stores
     * otherwise all of the sugar would just immediately go to the roots
     *
     * NOTE: the forces and therefore sugar transfered are worked out for all organs based on the amount of sugar there in the beginning
     * this could lead to a slight error, but should be corrected by the starch latter just have to imagine that all of the sugar
     * goes to the allocated organs simultaneously
     */

    conc_gradient concentration_gradient;
    // TODO: think of the units here!
    concentration_gradient.needles_to_phloem = ((sugar.needles+starch.needles)/needles_mass - (sugar.phloem+starch.phloem)/7.410537931)/parameters.resistance_needles_to_phloem;
    concentration_gradient.phloem_to_roots = ((sugar.phloem+starch.phloem)/7.410537931 - (sugar.roots+starch.roots)/2.8)/parameters.resistance_phloem_to_roots;
    concentration_gradient.phloem_to_xylem_sh = ((sugar.phloem+starch.phloem)/7.410537931 - (sugar.xylem_sh+starch.xylem_sh)/74.10537931)/parameters.resistance_phloem_to_xylem_sh;
    concentration_gradient.phloem_to_xylem_st = ((sugar.phloem+starch.phloem)/7.410537931 - (sugar.xylem_st+starch.xylem_st)/8.65862069)/parameters.resistance_phloem_to_xylem_st;
    concentration_gradient.roots_to_myco = parameters.mycorrhiza_threshold * (sugar.roots + starch.roots);

    /*
     * Balance calculations
     */

    double carbo_beginning = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
      starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + PF;
    double sugar_out_of_system = 0;

    /*
     * SUGAR TRANSFER WITH ALL PROCESSES BUT EMERGANCY
     */

    // # Rm.a maintenance respiration separated into organs
    sugar.needles = sugar.needles + PF -          // Last day sugar + daily photosynthesis
      resp.RmN * storage_term.needles -                          // maintenance respiration (altered by the carbon storage)
      (1 + common.Rg_N) * storage_term.needles * (pot_growth.needles + pot_growth.bud) -          // growth and growth respiration altered by the storage
      pot_growth.use + pot_growth.release -                                                           // growth sugar use and + release and to the rest of the organs
      concentration_gradient.needles_to_phloem*needles_mass +                            // transfer between organs
      (Kd.needles - Ks.needles) * parameters.carbon_sugar * 0.001 * needles_mass;   // + sperling processes with links to the needles growth process

    // coefficients are from mass ratio in starch and sugar 2015 xls
    sugar.phloem = sugar.phloem -
      0.082179938 * resp.RmS * storage_term.phloem -
      0.082179938 * (1 + common.Rg_S) * storage_term.phloem * (pot_growth.wall + pot_growth.height) +  // growth
      concentration_gradient.needles_to_phloem*needles_mass -                             // transfer between organs
      concentration_gradient.phloem_to_roots*7.410537931 -                                 // transfer between organs
      concentration_gradient.phloem_to_xylem_sh*7.410537931 -
      concentration_gradient.phloem_to_xylem_st*7.410537931 +
      (Kd.phloem - Ks.phloem) * parameters.carbon_sugar * 0.001 * 7.4;

    sugar.roots = sugar.roots -
      resp.RmR * storage_term.roots -                                              // maintenance respiration);
      (1 + common.Rg_R) * storage_term.roots * pot_growth.roots +             // growth
      concentration_gradient.phloem_to_roots*7.410537931 -        // transfer between organs
      concentration_gradient.roots_to_myco +                                         // transfer between organs, no multiplier as this is for mycorhiza and the model just takes the extra sugar
      (Kd.roots - Ks.roots) * parameters.carbon_sugar * 0.001 * 2.8;

    sugar.xylem_sh = sugar.xylem_sh -
      0.096020683 * resp.RmS * storage_term.xylem_sh -                                   // maintenance respiration
      0.096020683 * (1 + common.Rg_S) * storage_term.xylem_sh * (pot_growth.wall + pot_growth.height) +  // growth
      concentration_gradient.phloem_to_xylem_sh*7.410537931 +
      (Kd.xylem_sh - Ks.xylem_sh) * parameters.carbon_sugar * 0.001 * 8.65862069;

    sugar.xylem_st = sugar.xylem_st -
      0.821799379 * resp.RmS * storage_term.xylem_st -                                // maintenance respiration
      0.821799379 * (1 + common.Rg_S) * storage_term.xylem_st * (pot_growth.wall + pot_growth.height) +  // growth
      concentration_gradient.phloem_to_xylem_st*7.410537931 +
      (Kd.xylem_st - Ks.xylem_st) * parameters.carbon_sugar * 0.001 * 74.10537931;

    sugar.mycorrhiza = concentration_gradient.roots_to_myco;

    double respiration_growth = resp.RmN * storage_term.needles + (1 + common.Rg_N) * storage_term.needles * (pot_growth.needles + pot_growth.bud) +
      0.082179938 * resp.RmS * storage_term.phloem + 0.082179938 * (1 + common.Rg_S) * storage_term.phloem * (pot_growth.wall + pot_growth.height) +
      0.821799379 * resp.RmS * storage_term.xylem_st + 0.821799379 * (1 + common.Rg_S) * storage_term.xylem_st * (pot_growth.wall + pot_growth.height) +
      0.096020683 * resp.RmS * storage_term.xylem_sh + 0.096020683 * (1 + common.Rg_S) * storage_term.xylem_sh * (pot_growth.wall + pot_growth.height) +
      resp.RmR * storage_term.roots + (1 + common.Rg_R) * storage_term.roots * pot_growth.roots;

    /*
     * STARCH UPDATED SPERLING
     */

    // SPERLING MODEL

    starch.needles = starch.needles + (- Kd.needles + Ks.needles) * parameters.carbon_sugar * 0.001 * needles_mass; // Subtract starch degradation and add synthase to ST
    starch.phloem = starch.phloem + (- Kd.phloem + Ks.phloem) * parameters.carbon_sugar * 0.001 * 7.4; // Subtract starch degradation and add synthase to ST
    starch.roots = starch.roots + (- Kd.roots + Ks.roots) * parameters.carbon_sugar * 0.001 * 2.8; // Subtract starch degradation and add synthase to ST
    starch.xylem_sh = starch.xylem_sh + (- Kd.xylem_sh + Ks.xylem_sh) * parameters.carbon_sugar * 0.001 * 8.65862069; // Subtract starch degradation and add synthase to starch
    starch.xylem_st = starch.xylem_st + (- Kd.xylem_st + Ks.xylem_st) * parameters.carbon_sugar * 0.001 * 74.10537931; // Subtract starch degradation and add synthase to starch

    /*
     * STARCH AND SUGAR UPDATED EMERGANCY MODEL
     *
     * If sugar is below a certain value starch is released so the sugar doesn't go negative before starch
     *  This is a proxy for a starch metabolism system, which seems to be present under stress in literature
     *  but I can't find a mechanism for scots pine
     *  values are below the lowest recorded value
     */

    // Work out energy transfer

    carbo_tracker emergancy_transfer;
    emergancy_transfer.needles =  emergancy(sugar.needles, starch.needles, parameters.tau_emergancy_needles, parameters.lower_bound_needles);
    emergancy_transfer.phloem = emergancy(sugar.phloem, starch.phloem, parameters.tau_emergancy_phloem, parameters.lower_bound_phloem);
    emergancy_transfer.xylem_sh = emergancy(sugar.xylem_sh, starch.xylem_sh, parameters.tau_emergancy_xylem_sh, parameters.lower_bound_xylem_sh);
    emergancy_transfer.xylem_st = emergancy(sugar.xylem_st, starch.xylem_st, parameters.tau_emergancy_xylem_st, parameters.lower_bound_xylem_st);
    emergancy_transfer.roots = emergancy(sugar.roots, starch.roots, parameters.tau_emergancy_roots, parameters.lower_bound_roots);

    // sugar update

    sugar.needles = sugar.needles + emergancy_transfer.needles;
    sugar.phloem = sugar.phloem + emergancy_transfer.phloem;
    sugar.xylem_sh = sugar.xylem_sh + emergancy_transfer.xylem_sh;
    sugar.xylem_st = sugar.xylem_st + emergancy_transfer.xylem_st;
    sugar.roots = sugar.roots + emergancy_transfer.roots;

    // starch update

    starch.needles = starch.needles - emergancy_transfer.needles;
    starch.phloem = starch.phloem - emergancy_transfer.phloem;
    starch.xylem_sh = starch.xylem_sh - emergancy_transfer.xylem_sh;
    starch.xylem_st = starch.xylem_st - emergancy_transfer.xylem_st;
    starch.roots = starch.roots - emergancy_transfer.roots;

    /*
     * Mass check
     */

    sugar_out_of_system = respiration_growth + sugar.mycorrhiza + pot_growth.use - pot_growth.release;
    double carbo_ending = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
      starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + sugar_out_of_system;
    double difference = carbo_beginning - carbo_ending;
    if (difference > 0.00000000000001) { // 10^13
      std::cout << "On day " << day + 1 << " The carbohydrate balance compared to the beginning of the day " << difference << "\n";
    }

    /*
     * Storage check
     */

    if (sugar.needles <= 0 & starch.needles <= 0) {
      std::cerr << " Day " << day << " No Storage needles! Plant died" << "\n";
    }
    if (sugar.phloem <= 0 & starch.phloem <= 0) {
      std::cerr << " Day " << day << " No Storage phloem! Plant died" << "\n";
    }
    if (sugar.roots <= 0 & starch.roots <= 0) {
      std::cerr <<  " Day " << day << " No Storage roots! Plant died" << "\n";
    }
    if (sugar.xylem_sh <= 0 & starch.xylem_sh <= 0) {
      std::cerr << " Day " << day << " No Storage xylem shoot! Plant died" << "\n";
    }
    if (sugar.xylem_st <= 0 & starch.xylem_st <= 0) {
      std::cerr << " Day " << day << " No Storage xylem stem! Plant died" << "\n";
    }

    /*
     * SPERLING PARAMETER UPDATE FOR NEXT ITERATION
     */

    if (day == 0) {
      parameters_in.As = As_initiliser(parameters_in.Ad, temperature_equilibrium, parameters_in.Ad.B, parameters_in.As.B);
    }

    As_new.needles = (1-parameters.lambda_needles)*parameters_in.As.needles;
    As_new.phloem = (1-parameters.lambda_phloem)*parameters_in.As.phloem;
    As_new.roots = (1-parameters.lambda_roots)*parameters_in.As.roots;
    As_new.xylem_sh = (1-parameters.lambda_xylem_sh)*parameters_in.As.xylem_sh;
    As_new.xylem_st = (1-parameters.lambda_xylem_st)*parameters_in.As.xylem_st;

    Ad_new.needles = (1-parameters.lambda_needles)*parameters_in.Ad.needles;
    Ad_new.phloem = (1-parameters.lambda_phloem)*parameters_in.Ad.phloem;
    Ad_new.roots = (1-parameters.lambda_roots)*parameters_in.Ad.roots;
    Ad_new.xylem_sh = (1-parameters.lambda_xylem_sh)*parameters_in.Ad.xylem_sh;
    Ad_new.xylem_st = (1-parameters.lambda_xylem_st)*parameters_in.Ad.xylem_st;

    /*
     * Induce starch synthase if SC is high or degradation if it is low
     * These numbers are from september 2018
     * xylem not changed as no data to support it
     *
     * TODO: should I make these numbers automatic somehow?
     */

    if (sugar.needles>0.12) {
      As_new.needles=As_new.needles+parameters.delta_needles;
    }
    else if (sugar.needles<0.12 && starch.needles> 0) {
      Ad_new.needles=Ad_new.needles+parameters.delta_needles;
    }

    if  (sugar.phloem>0.28)
    {
      As_new.phloem=As_new.phloem+parameters.delta_phloem;
    }
    else if (sugar.phloem<0.28 && starch.phloem> 0)
    {
      Ad_new.phloem=Ad_new.phloem+parameters.delta_phloem;
    }

    if  (sugar.roots>0.09)
    {
      As_new.roots=As_new.roots+parameters.delta_roots;
    }
    else if (sugar.roots<0.09 && starch.roots> 0)
    {
      Ad_new.roots=Ad_new.roots+parameters.delta_roots;
    }

    if  (sugar.xylem_sh>0.049)
    {
      As_new.xylem_sh=As_new.xylem_sh+parameters.delta_xylem_sh;
    }
    else if (sugar.xylem_sh<0.049 && starch.xylem_sh > 0)
    {
      Ad_new.xylem_sh=Ad_new.xylem_sh+parameters.delta_xylem_sh;
    }

    if  (sugar.xylem_st>0.32)
    {
      As_new.xylem_st=As_new.xylem_st+parameters.delta_xylem_st;
    }
    else if (sugar.xylem_st<0.32 && starch.xylem_sh > 0)
    {
      Ad_new.xylem_st=Ad_new.xylem_st+parameters.delta_xylem_st;
    }

    /*
     * Bud burst trigger
     */
    // If it hasn't been bud burst yet then possible bud burst is calculated
    if (day < parameters.sB0-1) {
      if (sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots < parameters.SCb) {
        std::cout << "New Bud Burst set: " << day << "\n";
        sB0 = day;
      } else {
        sB0 = parameters.sB0;
      }
    }
  } else {
    // Sperling parameters
    As_new.needles = 0;
    As_new.phloem = 0;
    As_new.roots = 0;
    As_new.xylem_sh = 0;
    As_new.xylem_st = 0;
    Ad_new.needles = 0;
    Ad_new.phloem = 0;
    Ad_new.roots = 0;
    Ad_new.xylem_sh = 0;
    Ad_new.xylem_st = 0;
    sB0 = parameters.sB0;

    // Model
    double ak = 1 / (1 - 1/exp(parameters.alfa * (0.7430989 - parameters.Wala)));
    double storage, storage_term_Rm, sugar_all, starch_all, to_sugar, to_starch;
    double myco_allocation;
    if (day == 0) {
      sugar_all = sugar.needles = 0.4184208;
      starch_all = starch.needles = parameters.starch00;
      to_sugar = 0;
      to_starch = 0;
      storage = storage_term.respiration = 1;
    } else {
      storage = std::max(0.0 , std::min(1.0 , ak * (1.0 - 1.0 / exp(parameters.alfa * (sugar.needles + starch.needles - parameters.Wala)))));
      if (storage < 0.1) {
        storage_term.respiration = 0;
      } else {
        storage_term.respiration = 1;
      }

      if ((sH > parameters.sHc) & (sugar.needles + starch.needles > 0.07)) {
        myco_allocation = PF * 0.3;
      } else {
        myco_allocation = 0.0;
      }

      sugar_all = sugar.needles + PF - pot_growth.use + pot_growth.release - storage_term_Rm * resp.Rm_a -
        (1 + common.Rg_S) * storage * pot_growth.height -
        (1 + common.Rg_S) * storage * pot_growth.diameter -
        (1 + common.Rg_N) * storage * pot_growth.needles -
        (1 + common.Rg_R) * storage * pot_growth.roots -
        (1 + common.Rg_N) * storage * pot_growth.bud -
        myco_allocation;

      if (sugar_all < 0.41) {
        to_sugar = std::min(starch.needles, (0.41 - sugar_all) / parameters.tau_t);
        to_starch = 0;
      } else if (sugar_all > 0.41) {
        to_sugar = 0;
        to_starch = (sugar_all - 0.41) / parameters.tau_s;
      } else {
        to_sugar = 0;
        to_starch = 0;
      }
      starch_all = starch.needles + to_starch - to_sugar;
      sugar_all = sugar_all + to_sugar - to_starch;

      storage_term.needles = storage;
      storage_term.phloem = storage;
      storage_term.roots = storage;
      storage_term.xylem_sh = storage;
      storage_term.xylem_st = storage;
    }

    if (sugar.needles <= 0 & starch.needles <= 0) {
      std::cout << " No Storage! Plant died" << "\n";
      tree_alive = FALSE;
    }

    // I don't want to change the structure too much from the original sugar model
    // so the average sugar is reported as sugar needles
    sugar.needles = sugar_all;
    sugar.phloem = 0.0;
    sugar.roots = 0.0;
    sugar.xylem_sh = 0.0;
    sugar.xylem_st = 0.0;
    sugar.mycorrhiza = myco_allocation;

    starch.needles = starch_all;
    starch.phloem = 0.0;
    starch.roots = 0.0;
    starch.xylem_sh = 0.0;
    starch.xylem_st = 0.0;
  }

  /*
   * OUTPUT!
   */

  carbo_values_out previous_values_out;
  previous_values_out.storage_term = storage_term;
  previous_values_out.Ad = Ad_new;
  previous_values_out.As = As_new;
  previous_values_out.sB0 = sB0;
  previous_values_out.tree_alive = tree_alive;

  carbo_balance out;
  out.sugar = sugar;
  out.starch = starch;
  out.storage = storage_term;
  out.previous_values = previous_values_out;

  return out;
}
