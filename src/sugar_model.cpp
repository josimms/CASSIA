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
  if (!tree_alive) {
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

double nitrogen_storage(double sugar_mycorrhiza, std::string organ) {
  double out;
  /*
   Korhonen 2007, N mg g−1
   *
   Needles          4.9 ± 0.3
   Large branches   4.7 ± 0.1
   Small branches   4.4 ± 0.3
   Bark             3.7 ± 0.2
   Cones and seeds  2.1 ± 0.2
   *
   *
   Ding 2019
   *
   Root             TOOD
   *
   */

  if (organ == "needles") {
    out = std::max(0.0, 1.0/(1.0+exp(-3.0*(sugar_mycorrhiza - 4.0))));
  } else if (organ == "bud") {
    out = std::max(0.0, 1.0/(1.0+exp(-3.0*(sugar_mycorrhiza - 2.0))));
  } else if (organ == "wall") {
    out = std::max(0.0, 1.0/(1.0+exp(-3.0*(sugar_mycorrhiza - 3.0))));
  } else if (organ == "height") {
    out = std::max(0.0, 1.0/(1.0+exp(-3.0*(sugar_mycorrhiza - 3.0))));
  } else if (organ == "roots") {
    out = std::max(0.0, 1.0/(1.0+exp(-3.0*(sugar_mycorrhiza - 2.0)))); // TODO: what is the ratio here
  } else if (organ == "all") {
    out = std::max(0.0, 1.0/(1.0+exp(-3.0*(sugar_mycorrhiza - 4.0))));
  }

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

                          double nitrogen_balance,
                          bool nitrogen_change,
                          bool nitrogen_contrast,

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

  growth_out nitrogen_capacity;
  nitrogen_capacity.needles = nitrogen_storage(nitrogen_balance, "needles");
  nitrogen_capacity.bud = nitrogen_storage(nitrogen_balance, "bud");
  nitrogen_capacity.wall = nitrogen_storage(nitrogen_balance, "wall");
  nitrogen_capacity.height = nitrogen_storage(nitrogen_balance, "height");
  nitrogen_capacity.roots = nitrogen_storage(nitrogen_balance, "roots");

  double growth = 0;
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
       * Storage for respiration
       */

      carbo_tracker storage_term_resp;
      storage_term_resp.needles = respiration_check(storage_term.needles, parameters.lower_bound_needles);
      storage_term_resp.phloem = respiration_check(storage_term.phloem, parameters.lower_bound_phloem);
      storage_term_resp.xylem_sh = respiration_check(storage_term.xylem_sh, parameters.lower_bound_xylem_sh);
      storage_term_resp.xylem_st = respiration_check(storage_term.xylem_st, parameters.lower_bound_xylem_st);
      storage_term_resp.roots = respiration_check(storage_term.roots, parameters.lower_bound_roots);

      /*
       * The concentration gradient is represented by the difference in storage terms.
       *
       * This is then multiplied by the total storage in the organ and divided by a delay factor.
       *
       * The amount  that is transferred should then either be the sugar that is over the storage capacity of the organ
       */

      // NOTE various amount of organ mass is storage
      conc_gradient concentration_gradient;

      double needle_transfer = (sugar.needles + starch.needles) * (storage_term.needles - storage_term.phloem)/2;
      if (surplus_c) {
        needle_transfer = (sugar.needles + starch.needles) - 0.11 * needles_mass;
      }
      concentration_gradient.needles_to_phloem = std::max((sugar.needles + starch.needles) - 0.11 * needles_mass, needle_transfer);

      double phloem_capacity = (sugar.phloem + starch.phloem) - 0.11 * phloem_mass;

      double xylem_sh_transfer = (sugar.phloem + starch.phloem) * (storage_term.phloem - storage_term.xylem_sh)/2;
      if (surplus_c) {
        // NOTE: Phloem used here as the sugar is coming from the phloem
        xylem_sh_transfer = 0;
      }
      concentration_gradient.phloem_to_xylem_sh = std::max(xylem_sh_transfer, xylem_sh_respiration_share * phloem_capacity); // TODO; split somehow?

      double xylem_st_transfer = (sugar.phloem + starch.phloem) * (storage_term.phloem - storage_term.xylem_st)/2;
      if (surplus_c) {
        // NOTE: Phloem used here as the sugar is coming from the phloem
        xylem_st_transfer = 0;

      }
      concentration_gradient.phloem_to_xylem_st = std::max(xylem_st_transfer, xylem_st_respiration_share * phloem_capacity);

      double phloem_transfer = (sugar.phloem + starch.phloem) * (storage_term.phloem - storage_term.roots)/2;
      if (surplus_c) {
        phloem_transfer = phloem_respiration_share * phloem_capacity;
      }
      concentration_gradient.phloem_to_roots = std::max(phloem_respiration_share * phloem_capacity, phloem_transfer);

      double root_capacity = std::max(sugar.roots + starch.roots - 0.15 * root_mass, 0.0);
      double myco_transfer = parameters.mycorrhiza_threshold * (sugar.roots);
      if (surplus_c) {
        myco_transfer = root_capacity;
      } else if (nitrogen_contrast) {
        // Difference in total storage terms
        double nitrogen_capacity_average = (nitrogen_capacity.needles + nitrogen_capacity.bud + nitrogen_capacity.wall + nitrogen_capacity.height + nitrogen_capacity.roots)/5;
        double normalised_nitrogen_target_increase = std::max((storage_term.needles + storage_term.phloem + storage_term.roots + storage_term.xylem_sh + storage_term.xylem_st)/5 - nitrogen_capacity_average, 0.0);

        // Sugar needed to make this difference
        // TODO: replace this formula as this doesn't represent the uptake anymore
        double nitrogen_target_increase_sugar_investment = (1.0/4.0) * std::log((1.0/nitrogen_capacity_average - 1.0)/(1.0/(nitrogen_capacity_average+normalised_nitrogen_target_increase) - 1.0));

        // Total sugar and storage in the model
        double total_sugar = sugar.roots + sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st +
          starch.roots + starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st;
        // And the total sugar used considering only the sugar limitation - gives the sugar left over
        double total_use =             resp.RmR * storage_term_resp.roots     +                              (1.0 + common.Rg_R) * storage_term.roots    * pot_growth.roots +
                                       resp.RmN * storage_term_resp.needles   +                              (1.0 + common.Rg_N) * storage_term.needles  * (pot_growth.needles + pot_growth.bud) +
            phloem_respiration_share * resp.RmS * storage_term_resp.phloem    +   phloem_respiration_share * (1.0 + common.Rg_S) * storage_term.phloem   * (pot_growth.wall + pot_growth.height) +
          xylem_sh_respiration_share * resp.RmS * storage_term_resp.xylem_sh  + xylem_sh_respiration_share * (1.0 + common.Rg_S) * storage_term.xylem_sh * (pot_growth.wall + pot_growth.height) +
          xylem_st_respiration_share * resp.RmS * storage_term_resp.xylem_st  + xylem_st_respiration_share * (1.0 + common.Rg_S) * storage_term.xylem_st * (pot_growth.wall + pot_growth.height);

        double sugar_extra = std::max(total_sugar - total_use, 0.0);
        // Then considering the investment needed compared to the sugar in the roots. Note that the investment needed should always be positive
        myco_transfer = std::min(nitrogen_target_increase_sugar_investment, sugar_extra);
      }
      concentration_gradient.roots_to_myco = std::max(myco_transfer, root_capacity); // TODO; should this be the entire tree or the root?

      double xylem_sh_capacity = std::max((sugar.xylem_sh + starch.xylem_sh) - 0.002 * xylem_sh_mass, 0.0);
      double xylem_st_capacity = std::max((sugar.xylem_st + starch.xylem_st) - 0.05  * xylem_st_mass, 0.0);

      /*
       * Balance calculations
       */

      double carbo_beginning = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + PF;
      double sugar_out_of_system = 0;

      /*
       * Emzyme system
       */
      carbo_tracker enzyme;
      enzyme.needles  = std::max((Kd.needles - Ks.needles)   * parameters.carbon_sugar * 0.001 * needles_mass,  0.0);
      enzyme.phloem   = std::max((Kd.phloem - Ks.phloem)     * parameters.carbon_sugar * 0.001 * phloem_mass,   0.0);
      enzyme.xylem_sh = std::max((Kd.xylem_sh - Ks.xylem_sh) * parameters.carbon_sugar * 0.001 * xylem_sh_mass, 0.0);
      enzyme.xylem_st = std::max((Kd.xylem_st - Ks.xylem_st) * parameters.carbon_sugar * 0.001 * xylem_st_mass, 0.0);
      enzyme.roots    = std::max((Kd.roots - Ks.roots)       * parameters.carbon_sugar * 0.001 * root_mass,     0.0);

      /*
       * SUGAR TRANSFER WITH ALL PROCESSES BUT EMERGANCY
       */

      // # Rm.a maintenance respiration separated into organs
      sugar.needles = sugar.needles + PF -                                                          // Last day sugar + daily photosynthesis
        resp.RmN * storage_term_resp.needles -                                                           // maintenance respiration (altered by the carbon storage)
        (1 + common.Rg_N) * (std::min(storage_term.needles, nitrogen_capacity.needles) * pot_growth.needles + std::min(storage_term.needles, nitrogen_capacity.bud) * pot_growth.bud) -          // growth and growth respiration altered by the storage
        concentration_gradient.needles_to_phloem; // +                                                     // transfer between organs
        // enzyme.needles;                                                                                  // sperling processes with links to the needles growth process

      // coefficients are from mass ratio in starch and sugar 2015 xls
      sugar.phloem = sugar.phloem -
        phloem_respiration_share * resp.RmS * storage_term_resp.phloem -
        phloem_respiration_share * (1 + common.Rg_S) * (std::min(storage_term.phloem, nitrogen_capacity.wall) * pot_growth.wall + std::min(storage_term.phloem, nitrogen_capacity.height) * pot_growth.height) +  // growth
        concentration_gradient.needles_to_phloem -                                           // transfer between organs
        concentration_gradient.phloem_to_roots -                                             // transfer between organs
        concentration_gradient.phloem_to_xylem_sh -
        concentration_gradient.phloem_to_xylem_st +
        xylem_st_capacity +
        xylem_sh_capacity -
        pot_growth.use + pot_growth.release; // +                                                // growth sugar use and + release and to the rest of the organs
        // enzyme.phloem;

      sugar.roots = sugar.roots -
        resp.RmR * storage_term_resp.roots -                                                // maintenance respiration);
        (1 + common.Rg_R) * pot_growth.roots * std::min(storage_term.roots, nitrogen_capacity.roots) +                    // growth
        concentration_gradient.phloem_to_roots -                                           // transfer between organs
        concentration_gradient.roots_to_myco; // +                                             // transfer between organs, no multiplier as this is for mycorrhiza and the model just takes the extra sugar
        // enzyme.roots;

      sugar.xylem_sh = sugar.xylem_sh -
        xylem_sh_respiration_share * resp.RmS * storage_term_resp.xylem_sh -                                   // maintenance respiration
        xylem_sh_respiration_share * (1 + common.Rg_S) * (std::min(storage_term.xylem_sh, nitrogen_capacity.wall) * pot_growth.wall + std::min(storage_term.xylem_sh, nitrogen_capacity.height) * pot_growth.height) -               // growth
        xylem_sh_capacity +                                                                                  // Over the storage limit
        concentration_gradient.phloem_to_xylem_sh; // +
        // enzyme.xylem_sh;

      sugar.xylem_st = sugar.xylem_st -
        xylem_st_respiration_share * resp.RmS * storage_term_resp.xylem_st -                                // maintenance respiration
        xylem_st_respiration_share * (1 + common.Rg_S) * (std::min(storage_term.xylem_st, nitrogen_capacity.wall) * pot_growth.wall + std::min(storage_term.xylem_st, nitrogen_capacity.height) * pot_growth.height) -            // growth
        xylem_st_capacity +                                                                                // Over the storage limit
        concentration_gradient.phloem_to_xylem_st; // +
        // enzyme.xylem_st;

      /*
       * Respiration
       */

      respiration_growth =           (common.Rg_N) * (std::min(storage_term.needles, nitrogen_capacity.needles) * pot_growth.needles + std::min(storage_term.needles, nitrogen_capacity.bud)     * pot_growth.bud)    +
          phloem_respiration_share * (common.Rg_S) * (std::min(storage_term.phloem, nitrogen_capacity.wall)     * pot_growth.wall    + std::min(storage_term.phloem, nitrogen_capacity.height)   * pot_growth.height) +
        xylem_st_respiration_share * (common.Rg_S) * (std::min(storage_term.xylem_st, nitrogen_capacity.wall)   * pot_growth.wall    + std::min(storage_term.xylem_st, nitrogen_capacity.height) * pot_growth.height) +
        xylem_sh_respiration_share * (common.Rg_S) * (std::min(storage_term.xylem_sh, nitrogen_capacity.wall)   * pot_growth.wall    + std::min(storage_term.xylem_sh, nitrogen_capacity.height) * pot_growth.height) +
                                     (common.Rg_R) *  std::min(storage_term.roots, nitrogen_capacity.roots)     * pot_growth.roots;

      respiration_maintainence =        resp.RmN * storage_term_resp.needles  +
             phloem_respiration_share * resp.RmS * storage_term_resp.phloem   +
           xylem_sh_respiration_share * resp.RmS * storage_term_resp.xylem_sh +
           xylem_st_respiration_share * resp.RmS * storage_term_resp.xylem_st +
                                        resp.RmR * storage_term_resp.roots;

      growth =                        std::min(storage_term.needles, nitrogen_capacity.needles) * pot_growth.needles +
                                      std::min(storage_term.needles, nitrogen_capacity.bud)     * pot_growth.bud     +
          phloem_respiration_share * (std::min(storage_term.phloem, nitrogen_capacity.wall)     * pot_growth.wall    + std::min(storage_term.phloem, nitrogen_capacity.height)   * pot_growth.height) +
        xylem_st_respiration_share * (std::min(storage_term.xylem_st, nitrogen_capacity.wall)   * pot_growth.wall    + std::min(storage_term.xylem_st, nitrogen_capacity.height) * pot_growth.height) +
        xylem_sh_respiration_share * (std::min(storage_term.xylem_sh, nitrogen_capacity.wall)   * pot_growth.wall    + std::min(storage_term.xylem_sh, nitrogen_capacity.height) * pot_growth.height) +
                                      std::min(storage_term.roots, nitrogen_capacity.roots)     * pot_growth.roots;

      /*
       * Mycorrhiza
       */

      sugar.mycorrhiza = concentration_gradient.roots_to_myco;

      /*
       * Nitrogen
       */

      if (nitrogen_change) {
        nitrogen_capacity.needles = nitrogen_storage(nitrogen_balance/5.0, "needles");
        nitrogen_capacity.bud     = nitrogen_storage(nitrogen_balance/5.0, "bud");
        nitrogen_capacity.wall    = nitrogen_storage(nitrogen_balance/5.0, "wall");
        nitrogen_capacity.height  = nitrogen_storage(nitrogen_balance/5.0, "height");
        nitrogen_capacity.roots   = nitrogen_storage(nitrogen_balance/5.0, "roots");

        // TODO: actual uptake
        double nitrogen_uptake = sugar.mycorrhiza; // TODO: add an uptake function here!

        // C:N ratios are from the Korhonen 2013 paper
        nitrogen_balance = nitrogen_balance +
           nitrogen_uptake -
          (std::min(storage_term.needles, nitrogen_capacity.needles) * pot_growth.needles * 1.0/104.0     + std::min(storage_term.needles, nitrogen_capacity.bud)     * pot_growth.bud    * 1.0/221.0) -
          (std::min(storage_term.phloem, nitrogen_capacity.wall)     * pot_growth.wall    * 1.0/134.0     + std::min(storage_term.phloem, nitrogen_capacity.height)   * pot_growth.height * 1.0/134.0) -
          (std::min(storage_term.xylem_st, nitrogen_capacity.wall)   * pot_growth.wall    * 1.0/134.0     + std::min(storage_term.xylem_st, nitrogen_capacity.height) * pot_growth.height * 1.0/134.0) -
          (std::min(storage_term.xylem_sh, nitrogen_capacity.wall)   * pot_growth.wall    * 1.0/134.0     + std::min(storage_term.xylem_sh, nitrogen_capacity.height) * pot_growth.height * 1.0/134.0) -
           std::min(storage_term.roots, nitrogen_capacity.roots)     * pot_growth.roots   * 1.0/100.0; // TODO: real value for the roots
      }

      /*
       * STARCH UPDATED SPERLING
       */

      // SPERLING MODEL

      starch.needles  = starch.needles; //   - enzyme.needles;        // Subtract starch degradation and add synthase to starch
      starch.phloem   = starch.phloem; //    - enzyme.phloem;         // Subtract starch degradation and add synthase to starch
      starch.roots    = starch.roots; //     - enzyme.roots;          // Subtract starch degradation and add synthase to starch
      starch.xylem_sh = starch.xylem_sh; //  - enzyme.xylem_sh;       // Subtract starch degradation and add synthase to starch
      starch.xylem_st = starch.xylem_st; //  - enzyme.xylem_st;       // Subtract starch degradation and add synthase to starch

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
      emergancy_transfer.needles  = emergancy(sugar.needles,  starch.needles,  parameters.tau_emergancy_needles,  parameters.lower_bound_needles);
      emergancy_transfer.phloem   = emergancy(sugar.phloem,   starch.phloem,   parameters.tau_emergancy_phloem,   parameters.lower_bound_phloem);
      emergancy_transfer.xylem_sh = emergancy(sugar.xylem_sh, starch.xylem_sh, parameters.tau_emergancy_xylem_sh, parameters.lower_bound_xylem_sh);
      emergancy_transfer.xylem_st = emergancy(sugar.xylem_st, starch.xylem_st, parameters.tau_emergancy_xylem_st, parameters.lower_bound_xylem_st);
      emergancy_transfer.roots    = emergancy(sugar.roots,    starch.roots,    parameters.tau_emergancy_roots,    parameters.lower_bound_roots);

      // Sugar update

      sugar.needles   = sugar.needles   + emergancy_transfer.needles;
      sugar.phloem    = sugar.phloem    + emergancy_transfer.phloem;
      sugar.xylem_sh  = sugar.xylem_sh  + emergancy_transfer.xylem_sh;
      sugar.xylem_st  = sugar.xylem_st  + emergancy_transfer.xylem_st;
      sugar.roots     = sugar.roots     + emergancy_transfer.roots;

      // Starch update

      starch.needles  = starch.needles   - emergancy_transfer.needles;
      starch.phloem   = starch.phloem    - emergancy_transfer.phloem;
      starch.xylem_sh = starch.xylem_sh  - emergancy_transfer.xylem_sh;
      starch.xylem_st = starch.xylem_st  - emergancy_transfer.xylem_st;
      starch.roots    = starch.roots     - emergancy_transfer.roots;

      /*
       * Mass check
       */

      sugar_out_of_system = growth + respiration_growth + respiration_maintainence + sugar.mycorrhiza + pot_growth.use - pot_growth.release;
      double carbo_ending = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + sugar_out_of_system;
      double difference = carbo_beginning - carbo_ending;
      if (difference > 0.00000000000001) { // 10^13
        std::cout << "On day " << day + 1 << " The carbohydrate balance compared to the beginning of the day " << difference << "\n";
      }

      /*
       * Storage check
       */

      if ((starch.needles == 0) && (starch.phloem == 0) && (starch.xylem_sh == 0) && (starch.xylem_st == 0) && (starch.roots == 0) &&
          (sugar.needles == 0) && (sugar.phloem == 0) && (sugar.xylem_sh == 0) && (sugar.xylem_st == 0) && (sugar.roots == 0)) {
        std::cerr << " Day " << day << " No total storage - plant died!" << "\n";
        // tree_alive = false;
      }

      /*
       * SPERLING PARAMETER UPDATE FOR NEXT ITERATION
       */

      if (day == 0) {
        parameters_in.As = As_initiliser(parameters_in.Ad, temperature_equilibrium, parameters_in.Ad.B, parameters_in.As.B);
      }

      As_new.needles  = (1-parameters.lambda_needles)  * parameters_in.As.needles;
      As_new.phloem   = (1-parameters.lambda_phloem)   * parameters_in.As.phloem;
      As_new.roots    = (1-parameters.lambda_roots)    * parameters_in.As.roots;
      As_new.xylem_sh = (1-parameters.lambda_xylem_sh) * parameters_in.As.xylem_sh;
      As_new.xylem_st = (1-parameters.lambda_xylem_st) * parameters_in.As.xylem_st;

      Ad_new.needles  = (1-parameters.lambda_needles)  * parameters_in.Ad.needles;
      Ad_new.phloem   = (1-parameters.lambda_phloem)   * parameters_in.Ad.phloem;
      Ad_new.roots    = (1-parameters.lambda_roots)    * parameters_in.Ad.roots;
      Ad_new.xylem_sh = (1-parameters.lambda_xylem_sh) * parameters_in.Ad.xylem_sh;
      Ad_new.xylem_st = (1-parameters.lambda_xylem_st) * parameters_in.Ad.xylem_st;

      /*
       * Induce starch synthase if SC is high or degradation if it is low
       * These numbers are from september 2018
       * xylem not changed as no data to support it
       *
       * TODO: should I make these numbers automatic somehow?
       */

      if (sugar.needles > parameters.sugar_needles00) {
        As_new.needles=As_new.needles+parameters.delta_needles;
      }
      else if (sugar.needles < parameters.sugar_needles00 && starch.needles > 0) {
        Ad_new.needles=Ad_new.needles+parameters.delta_needles;
      }

      if  (sugar.phloem > parameters.sugar_phloem00)
      {
        As_new.phloem=As_new.phloem+parameters.delta_phloem;
      }
      else if (sugar.phloem < parameters.sugar_phloem00 && starch.phloem > 0)
      {
        Ad_new.phloem=Ad_new.phloem+parameters.delta_phloem;
      }

      if (sugar.roots > parameters.sugar_roots00)
      {
        As_new.roots=As_new.roots+parameters.delta_roots;
      }
      else if (sugar.roots < parameters.sugar_roots00 && starch.roots > 0)
      {
        Ad_new.roots=Ad_new.roots+parameters.delta_roots;
      }

      if  (sugar.xylem_sh > parameters.sugar_xylem_sh00)
      {
        As_new.xylem_sh=As_new.xylem_sh+parameters.delta_xylem_sh;
      }
      else if (sugar.xylem_sh < parameters.sugar_xylem_sh00 && starch.xylem_sh > 0)
      {
        Ad_new.xylem_sh=Ad_new.xylem_sh+parameters.delta_xylem_sh;
      }

      if  (sugar.xylem_st > parameters.sugar_xylem_st00)
      {
        As_new.xylem_st=As_new.xylem_st+parameters.delta_xylem_st;
      }
      else if (sugar.xylem_st < parameters.sugar_xylem_st00 && starch.xylem_st > 0)
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
    double storage, storage_term_Rm, sugar_all, starch_all, to_sugar, to_starch, nitrogen_capacity;
    double myco_allocation;
    if (day == 0) {
      sugar_all = parameters.sugar0;
      starch_all = parameters.starch0;
      to_sugar = 0;
      to_starch = 0;
      storage = storage_term.respiration = 1;
      nitrogen_capacity = 1;
    } else {
      storage = std::max(0.0 , std::min(1.0 , ak * (1.0 - 1.0 / exp(parameters.alfa * (sugar.needles + starch.needles - parameters.Wala)))));
      nitrogen_capacity = nitrogen_storage(nitrogen_balance/5.0, "all");
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

    if ((sugar.needles <= 0) && (starch.needles <= 0)) {
      std::cout << "Day: " << day << " No Storage! Plant died" << "\n";
      // tree_alive = false;
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
    respiration_maintainence = storage_term.respiration * resp.Rm_a;
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
  out.nitrogen_capacity = nitrogen_capacity;
  out.resp_growth = respiration_growth;
  out.resp_main = respiration_maintainence;
  out.nitrogen_balance = nitrogen_balance;
  out.previous_values = previous_values_out;

  return out;
}
