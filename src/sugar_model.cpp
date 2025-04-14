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

// TODO: make this make sense for each organ
double storage_update_organs(double storage_capacity, double sugar, double starch, bool tree_alive) {
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
    out = std::max(0.0, 1/(1+exp(-2*(sugar + starch - storage_capacity/2))));
    if (out > 1) {
      out = 1; // Corrected here as this is calculated at the beginning, but in the model the extra sugar would be transferred at the end of the iteration
    }
  }
  return out;
}

double nitrogen_storage(double sugar_mycorrhiza) {
  double out;
  out = std::max(0.0, 1/(1+exp(-4*(sugar_mycorrhiza - 0.01))));
  return(out);
}

carbo_tracker As_initiliser(carbo_tracker Ad, carbo_tracker equilibrium_temperature, double Bd, double Bs)
{
  carbo_tracker As;
  As.needles = Ad.needles * std::exp(equilibrium_temperature.needles * (Bd - Bs));
  As.phloem = Ad.phloem * std::exp(equilibrium_temperature.phloem * (Bd - Bs));
  As.xylem_sh = Ad.xylem_sh * std::exp(equilibrium_temperature.xylem_sh * (Bd - Bs));
  As.xylem_st = Ad.xylem_st * std::exp(equilibrium_temperature.xylem_st * (Bd - Bs));
  As.roots = Ad.roots * std::exp(equilibrium_temperature.roots * (Bd - Bs));

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

double respiration_check(double storage_term, double lower_bound) {
  double out = storage_term;
  if (storage_term < lower_bound) {
    out = 0;
  }
  return(out);
}

/*
 * Some of the parameters are site dependent, could define the site side in the main function code so I just import the right variables into this code
 *
 * needles mass per year, should define the year before this function!
 */

carbo_balance sugar_model(int year,
                          int day,
                          double TAir,
                          double PF,

                          CASSIA_common common,
                          CASSIA_parameters parameters,

                          double D00,
                          double sH,
                          respiration_out resp,

                          double nitrogen_capacity,
                          bool nitrogen_change,

                          bool sperling_sugar_model,
                          bool tree_alive,
                          bool storage_grows,
                          bool surplus_c,
                          double needles_mass, // Repola
                          double root_mass,
                          carbo_tracker temperature_equilibrium, // Calculated in the main function

                          growth_out pot_growth,

                          carbo_tracker sugar,
                          carbo_tracker starch,

                          carbo_values_out parameters_in) {

  // TODO: phloem mass is hard to get within the model, maybe a product of diameter and height
  double phloem_mass = 7.410537931;
  double xylem_sh_mass = 74.10537931;
  double xylem_st_mass = 8.65862069;

  double phloem_respiration_share = phloem_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);
  double xylem_sh_respiration_share = xylem_sh_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);
  double xylem_st_respiration_share = xylem_st_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);

  carbo_tracker storage_term;
  // 7.5% of mass should be for storage, von Arx 2017 max axis
  storage_term.needles = storage_update_organs(0.11*needles_mass, sugar.needles, starch.needles, tree_alive);
  storage_term.phloem = storage_update_organs(0.11*phloem_mass, sugar.phloem, starch.phloem, tree_alive);
  storage_term.xylem_sh = storage_update_organs(0.002*xylem_sh_mass, sugar.xylem_sh, starch.xylem_sh, tree_alive);
  storage_term.xylem_st = storage_update_organs(0.05*xylem_st_mass, sugar.xylem_st, starch.xylem_st, tree_alive);
  storage_term.roots = storage_update_organs(0.15*root_mass, sugar.roots, starch.roots, tree_alive);

  double respiration_growth = 0;
  double respiration_maintainence = 0;

  double sB0;
  carbo_tracker As_new, Ad_new;

  if (sperling_sugar_model) {
    if (day == 0) {
      sugar.needles = parameters.sugar_needles0;
      sugar.phloem = parameters.sugar_phloem0;
      sugar.roots = parameters.sugar_roots0;
      sugar.xylem_sh = parameters.sugar_xylem_sh0;
      sugar.xylem_st = parameters.sugar_xylem_st0;
      sugar.mycorrhiza = parameters.mycorrhiza_threshold * sugar.roots;

      starch.needles = parameters.starch_needles0;
      starch.phloem = parameters.starch_phloem0;
      starch.roots = parameters.starch_roots0;
      starch.xylem_sh = parameters.starch_xylem_sh0;
      starch.xylem_st = parameters.starch_xylem_st0;
      starch.mycorrhiza = parameters.mycorrhiza_threshold * starch.roots;
    } else {
      /*
       * PARAMETER UPDATES!
       */

      carbo_tracker Kd = carbo_tracker_init();
      carbo_tracker Ks = carbo_tracker_init();

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
       * The concentration gradient is represented by the difference in storage terms.
       *
       * This is then multiplied by the total storage in the organ and divided by a delay factor.
       *
       * The amount  that is transferred should then either be the sugar that is over the storage capacity of the organ
       */

      // NOTE 7.5% of organ mass is storage
      conc_gradient concentration_gradient;

      double needle_transfer = (sugar.needles + starch.needles) * (storage_term.needles - storage_term.phloem)/2;
      concentration_gradient.needles_to_phloem = std::max((sugar.needles + starch.needles) - 0.11 * needles_mass, needle_transfer);

      double phloem_transfer = (sugar.phloem + starch.phloem) * (storage_term.phloem - storage_term.roots)/2;
      concentration_gradient.phloem_to_roots = std::max((sugar.phloem + starch.phloem) - 0.11 * phloem_mass, phloem_transfer);

      double xylem_sh_transfer = (sugar.phloem + starch.phloem) * (storage_term.phloem - storage_term.xylem_sh)/2;
      concentration_gradient.phloem_to_xylem_sh = xylem_sh_transfer;

      double xylem_st_transfer = (sugar.phloem + starch.phloem) * (storage_term.phloem - storage_term.xylem_st)/2;
      concentration_gradient.phloem_to_xylem_st = xylem_st_transfer;

      double root_capacity = std::max(sugar.roots + starch.roots - 0.15 * root_mass, 0.0);
      double myco_transfer = parameters.mycorrhiza_threshold * (sugar.roots + starch.roots);
      if (surplus_c) {
        myco_transfer = std::max(sugar.roots - resp.RmR * storage_term.roots - (1 + common.Rg_R) * pot_growth.roots, 0.0);
      }
      /*
      else if (nitrogen_contrast) {


        myco_transfer = ;
      }
       */
      concentration_gradient.roots_to_myco = std::max(myco_transfer, root_capacity);

      double xylem_sh_capacity = std::max((sugar.xylem_sh + starch.xylem_sh) - 0.002 * xylem_sh_mass, 0.0);
      double xylem_st_capacity = std::max((sugar.xylem_st + starch.xylem_st) - 0.05 * xylem_st_mass, 0.0);

      /*
       * Balance calculations
       */

      double carbo_beginning = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + PF;
      double sugar_out_of_system = 0;

      /*
       * Storage for respiration
       */

      carbo_tracker storage_term_resp;
      storage_term_resp.needles = respiration_check(storage_term.needles, parameters.lower_bound_needles);
      storage_term_resp.phloem = respiration_check(storage_term.phloem, parameters.lower_bound_phloem);
      storage_term_resp.xylem_sh = respiration_check(storage_term.xylem_sh, parameters.lower_bound_xylem_sh);
      storage_term_resp.xylem_st = respiration_check(storage_term.xylem_st, parameters.lower_bound_xylem_st);
      storage_term_resp.roots = respiration_check(storage_term.roots, parameters.lower_bound_roots);

      /*
       * SUGAR TRANSFER WITH ALL PROCESSES BUT EMERGANCY
       */

      // # Rm.a maintenance respiration separated into organs
      sugar.needles = sugar.needles + PF -                                                          // Last day sugar + daily photosynthesis
        resp.RmN * storage_term_resp.needles -                                                           // maintenance respiration (altered by the carbon storage)
        (1 + common.Rg_N) * std::min(storage_term.needles, nitrogen_capacity) * (pot_growth.needles + pot_growth.bud) -          // growth and growth respiration altered by the storage
        concentration_gradient.needles_to_phloem +                                                  // transfer between organs
        (Kd.needles - Ks.needles) * parameters.carbon_sugar * 0.001 * needles_mass;                 // + sperling processes with links to the needles growth process

      sugar.needles = std::max(sugar.needles, 0.0);

      // coefficients are from mass ratio in starch and sugar 2015 xls
      sugar.phloem = sugar.phloem -
        phloem_respiration_share * resp.RmS * storage_term_resp.phloem -
        phloem_respiration_share * (1 + common.Rg_S) * std::min(storage_term.phloem, nitrogen_capacity) * (pot_growth.wall + pot_growth.height) +  // growth
        concentration_gradient.needles_to_phloem -                                         // transfer between organs
        concentration_gradient.phloem_to_roots -                                           // transfer between organs
        concentration_gradient.phloem_to_xylem_sh -
        concentration_gradient.phloem_to_xylem_st +
        xylem_st_capacity +
        xylem_sh_capacity -
        pot_growth.use + pot_growth.release +                                              // growth sugar use and + release and to the rest of the organs
        (Kd.phloem - Ks.phloem) * parameters.carbon_sugar * 0.001 * phloem_mass;

      sugar.phloem = std::max(sugar.phloem, 0.0);

      sugar.roots = sugar.roots -
        resp.RmR * storage_term_resp.roots -                                                // maintenance respiration);
        (1 + common.Rg_R) * pot_growth.roots * std::min(storage_term.roots, nitrogen_capacity) +                    // growth
        concentration_gradient.phloem_to_roots -                                       // transfer between organs
        concentration_gradient.roots_to_myco +                                         // transfer between organs, no multiplier as this is for mycorrhiza and the model just takes the extra sugar
        (Kd.roots - Ks.roots) * parameters.carbon_sugar * 0.001 * root_mass;

      sugar.roots = std::max(sugar.roots, 0.0);

      sugar.xylem_sh = sugar.xylem_sh -
        xylem_sh_capacity -                                                                               // Over the storage limit
        xylem_sh_respiration_share * resp.RmS * storage_term_resp.xylem_sh -                                   // maintenance respiration
        xylem_sh_respiration_share * (1 + common.Rg_S) * std::min(storage_term.xylem_sh, nitrogen_capacity) * (pot_growth.wall + pot_growth.height) +               // growth
        concentration_gradient.phloem_to_xylem_sh +
        (Kd.xylem_sh - Ks.xylem_sh) * parameters.carbon_sugar * 0.001 * xylem_sh_mass;

      sugar.xylem_sh = std::max(sugar.xylem_sh, 0.0);

      sugar.xylem_st = sugar.xylem_st -
        xylem_st_capacity -                                                                            // Over the storage limit
        xylem_st_respiration_share * resp.RmS * storage_term_resp.xylem_st -                                // maintenance respiration
        xylem_st_respiration_share * (1 + common.Rg_S) * std::min(storage_term.xylem_st, nitrogen_capacity) * (pot_growth.wall + pot_growth.height) +            // growth
        concentration_gradient.phloem_to_xylem_st +
        (Kd.xylem_st - Ks.xylem_st) * parameters.carbon_sugar * 0.001 * xylem_st_mass;

      sugar.xylem_st = std::max(sugar.xylem_st, 0.0);

      /*
       * Respiration
       */

      respiration_growth = (common.Rg_N) * std::min(storage_term.needles, nitrogen_capacity) * (pot_growth.needles + pot_growth.bud) +
          phloem_respiration_share * (common.Rg_S) * std::min(storage_term.phloem, nitrogen_capacity) * (pot_growth.wall + pot_growth.height) +
        xylem_st_respiration_share * (common.Rg_S) * std::min(storage_term.xylem_st, nitrogen_capacity) * (pot_growth.wall + pot_growth.height) +
        xylem_sh_respiration_share * (common.Rg_S) * std::min(storage_term.xylem_sh, nitrogen_capacity) * (pot_growth.wall + pot_growth.height) +
                                     (common.Rg_R) * std::min(storage_term.roots, nitrogen_capacity) * pot_growth.roots;

      respiration_maintainence = resp.RmN * storage_term_resp.needles * nitrogen_capacity +
             phloem_respiration_share * resp.RmS * storage_term_resp.phloem * nitrogen_capacity +
           xylem_sh_respiration_share * resp.RmS * storage_term_resp.xylem_sh * nitrogen_capacity +
           xylem_st_respiration_share * resp.RmS * storage_term_resp.xylem_st * nitrogen_capacity +
                                        resp.RmR * storage_term_resp.roots * nitrogen_capacity;

      /*
       * Mycorrhiza
       */

      sugar.mycorrhiza = concentration_gradient.roots_to_myco;

      if (nitrogen_change) {
        nitrogen_capacity = nitrogen_storage(sugar.mycorrhiza);
      }

      /*
       * STARCH UPDATED SPERLING
       */

      // SPERLING MODEL

      starch.needles = starch.needles + (- Kd.needles + Ks.needles) * parameters.carbon_sugar * 0.001 * needles_mass;      // Subtract starch degradation and add synthase to ST
      starch.needles = std::max(starch.needles, 0.0);
      starch.phloem = starch.phloem + (- Kd.phloem + Ks.phloem) * parameters.carbon_sugar * 0.001 * phloem_mass;           // Subtract starch degradation and add synthase to ST
      starch.phloem = std::max(starch.phloem, 0.0);
      starch.roots = starch.roots + (- Kd.roots + Ks.roots) * parameters.carbon_sugar * 0.001 * root_mass;                 // Subtract starch degradation and add synthase to ST
      starch.roots = std::max(starch.roots, 0.0);
      starch.xylem_sh = starch.xylem_sh + (- Kd.xylem_sh + Ks.xylem_sh) * parameters.carbon_sugar * 0.001 * xylem_sh_mass; // Subtract starch degradation and add synthase to starch
      starch.xylem_sh = std::max(starch.xylem_sh, 0.0);
      starch.xylem_st = starch.xylem_st + (- Kd.xylem_st + Ks.xylem_st) * parameters.carbon_sugar * 0.001 * xylem_st_mass; // Subtract starch degradation and add synthase to starch
      starch.xylem_st = std::max(starch.xylem_st, 0.0);

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

      // Sugar update

      sugar.needles = sugar.needles + emergancy_transfer.needles;
      sugar.phloem = sugar.phloem + emergancy_transfer.phloem;
      sugar.xylem_sh = sugar.xylem_sh + emergancy_transfer.xylem_sh;
      sugar.xylem_st = sugar.xylem_st + emergancy_transfer.xylem_st;
      sugar.roots = sugar.roots + emergancy_transfer.roots;

      // Starch update

      starch.needles = starch.needles - emergancy_transfer.needles;
      starch.phloem = starch.phloem - emergancy_transfer.phloem;
      starch.xylem_sh = starch.xylem_sh - emergancy_transfer.xylem_sh;
      starch.xylem_st = starch.xylem_st - emergancy_transfer.xylem_st;
      starch.roots = starch.roots - emergancy_transfer.roots;

      /*
       * Mass check
       */

      sugar_out_of_system = respiration_growth + respiration_maintainence + sugar.mycorrhiza + pot_growth.use - pot_growth.release;
      double carbo_ending = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + sugar_out_of_system;
      double difference = carbo_beginning - carbo_ending;
      if (difference > 0.00000000000001) { // 10^13
        std::cout << "On day " << day + 1 << " The carbohydrate balance compared to the beginning of the day " << difference << "\n";
      }

      /*
       * Storage check
       */

      if ((starch.needles == 0) && (starch.phloem == 0) && (starch.xylem_sh == 0) && (starch.xylem_st == 0) && (starch.roots == 0)) {
        std::cerr << " Day " << day << " No total storage - plant died!" << "\n";
        tree_alive = false;
      }
      if ((sugar.needles == 0) && (sugar.phloem == 0) && (sugar.xylem_sh == 0) && (sugar.xylem_st == 0) && (sugar.roots == 0)) {
        std::cerr << " Day " << day << " No total sugar - plant died!" << "\n";
        tree_alive = false;
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

      if (sugar.needles > parameters.sugar_needles0) {
        As_new.needles=As_new.needles+parameters.delta_needles;
      }
      else if (sugar.needles < parameters.sugar_needles0 && starch.needles> 0) {
        Ad_new.needles=Ad_new.needles+parameters.delta_needles;
      }

      if  (sugar.phloem > parameters.sugar_phloem0)
      {
        As_new.phloem=As_new.phloem+parameters.delta_phloem;
      }
      else if (sugar.phloem < parameters.sugar_phloem0 && starch.phloem> 0)
      {
        Ad_new.phloem=Ad_new.phloem+parameters.delta_phloem;
      }

      if (sugar.roots > parameters.sugar_roots0)
      {
        As_new.roots=As_new.roots+parameters.delta_roots;
      }
      else if (sugar.roots < parameters.sugar_roots0 && starch.roots> 0)
      {
        Ad_new.roots=Ad_new.roots+parameters.delta_roots;
      }

      if  (sugar.xylem_sh > parameters.sugar_xylem_sh0)
      {
        As_new.xylem_sh=As_new.xylem_sh+parameters.delta_xylem_sh;
      }
      else if (sugar.xylem_sh < parameters.sugar_xylem_sh0 && starch.xylem_sh > 0)
      {
        Ad_new.xylem_sh=Ad_new.xylem_sh+parameters.delta_xylem_sh;
      }

      if  (sugar.xylem_st > parameters.sugar_xylem_st0)
      {
        As_new.xylem_st=As_new.xylem_st+parameters.delta_xylem_st;
      }
      else if (sugar.xylem_st < parameters.sugar_xylem_st0 && starch.xylem_sh > 0)
      {
        Ad_new.xylem_st=Ad_new.xylem_st+parameters.delta_xylem_st;
      }

      /*
       * Bud burst trigger
       */
      // If it hasn't been bud burst yet then possible bud burst is calculated
      /*
      if (day < parameters.sB0-1) {
        if (sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots < parameters.SCb) {
          std::cout << "New Bud Burst set: " << day << "\n";
          sB0 = day;
        } else {
          sB0 = parameters.sB0;
        }
      }
       */
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
    double ak = 1 / (1 - 1/exp(parameters.alfa * (parameters.sugar00 + parameters.starch00 - parameters.Wala)));
    double storage, storage_term_Rm, sugar_all, starch_all, to_sugar, to_starch;
    double myco_allocation;
    if (day == 0) {
      sugar_all = parameters.sugar0;
      starch_all = parameters.starch0;
      to_sugar = 0;
      to_starch = 0;
      storage = storage_term.respiration = 1;
    } else {
      storage = std::max(0.0 , std::min(1.0 , ak * (1.0 - 1.0 / exp(parameters.alfa * (sugar.needles + starch.needles - parameters.Wala)))));
      if ((sugar.needles + starch.needles) < 0.1) {
        storage_term.respiration = 0.0;
      } else {
        storage_term.respiration = 1.0;
      }

      if ((sH > parameters.sHc) & ((sugar.needles + starch.needles) > 0.07)) {
        myco_allocation = PF * 0.3;
      } else {
        myco_allocation = 0.0;
      }

      sugar_all = sugar.needles + PF - pot_growth.use + pot_growth.release - std::min(storage_term.respiration, nitrogen_capacity) * resp.Rm_a -
        (1 + common.Rg_S) * std::min(storage, nitrogen_capacity) * pot_growth.height -
        (1 + common.Rg_S) * std::min(storage, nitrogen_capacity) * pot_growth.diameter -
        (1 + common.Rg_N) * std::min(storage, nitrogen_capacity) * pot_growth.needles -
        (1 + common.Rg_R) * std::min(storage, nitrogen_capacity) * pot_growth.roots -
        (1 + common.Rg_N) * std::min(storage, nitrogen_capacity) * pot_growth.bud -
        myco_allocation;

      if (sugar_all < parameters.sugar00) {
        to_sugar = std::min(starch.needles, (parameters.sugar00 - sugar_all) / 2.0);
        to_starch = 0;
      } else if (sugar_all > parameters.sugar00) {
        to_sugar = 0;
        to_starch = (sugar_all - parameters.sugar00) / 2.0;
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

    if ((sugar.needles <= 0) & (starch.needles <= 0)) {
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

    respiration_growth = (common.Rg_S) * std::min(storage, nitrogen_capacity) * pot_growth.height -
      (common.Rg_S) * std::min(storage, nitrogen_capacity) * pot_growth.diameter -
      (common.Rg_N) * std::min(storage, nitrogen_capacity) * pot_growth.needles -
      (common.Rg_R) * std::min(storage, nitrogen_capacity) * pot_growth.roots -
      (common.Rg_N) * std::min(storage, nitrogen_capacity) * pot_growth.bud;
    respiration_maintainence = std::min(storage_term.respiration, nitrogen_capacity) * resp.Rm_a;
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
  out.resp_growth = respiration_growth;
  out.resp_main = respiration_maintainence;
  out.nitrogen_capacity = nitrogen_capacity;
  out.previous_values = previous_values_out;

  return out;
}
