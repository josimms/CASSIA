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

double storage_update_organs(double growth, double sugar, double lower_limit, bool tree_alive) {
  if (!tree_alive || std::isnan(sugar)) {
    return 0.0;
  }
  double out;

  // Sugar available beyond the minimum survival level
  double available = sugar - lower_limit;

  if (available >= 0) {
    // Avoid division by zero in case growth is very small
    double ratio = available / (growth + 1e-6);

    // Smooth sigmoid scaling — center at ratio = 1
    // double scale = 1.0 / (1.0 + std::exp(-4.0 * (ratio - 1.0)));

    // Final result: limited to 1, capped so usage won't exceed available sugar
    out = std::min(1.0, ratio);
  } else {
    out = 0.0;
  }

  return out;
}


double sugar_demand(double input) {
  double out = 1.0/(1.0 + std::exp(-1*input));
  return out;
}

double safe_to_starch_transfer(
    double sugar,
    double winter,
    double starch,
    double Ad0,
    double delta,
    double space_left  // Available starch storage capacity
) {
  // Sugar surplus after winter needs
  double S = sugar - winter;

  if (S <= 0.0) {
    // Sugar is low, remobilize from starch

    double sd_sugar = sugar_demand(S);  // likely negative or near zero
    double sd_starch = sugar_demand(starch);

    double proposed_flux = Ad0 * 0.0 - delta * (1.0 - sd_sugar) * std::abs(S) * sd_starch;

    // Cap by available starch (negative flux)
    proposed_flux = std::max(proposed_flux, -starch);

    return proposed_flux;
  }

  // Sugar surplus available to convert to starch
  double sd_sugar  = sugar_demand(S);
  if (sd_sugar > 1.0 || sd_sugar < 0.0) {
    std::cout << "sd_sugar is < 0.0 or > 1.0\n";
  }
  double sd_starch = sugar_demand(starch);
  if (sd_starch > 1.0 || sd_starch < 0.0) {
    std::cout << "sd_starch is < 0.0 or > 1.0\n";
  }

  // Calculate proposed flux sugar → starch (positive flux)
  double proposed_flux = Ad0 * sd_sugar * S - delta * (1.0 - sd_sugar) * S * sd_starch;

  if (proposed_flux > 0.0) {
    // Cap by sugar availability AND starch storage capacity
    proposed_flux = std::min(proposed_flux, sugar);
    proposed_flux = std::min(proposed_flux, space_left);
  } else {
    // Cap by starch availability (for starch → sugar)
    proposed_flux = std::max(proposed_flux, -starch);
  }

  return proposed_flux;
}


double nitrogen_storage(
    double nitrogen_balance,
    std::string organ
) {

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


  double storage_capcity = 0.0;
  if (organ == "needles") {
    storage_capcity = 4.9;
  } else if (organ == "bud") {
    storage_capcity = 2.1;
  } else if (organ == "wall") {
    storage_capcity = 3.7;
  } else if (organ == "height") {
    storage_capcity = 4.7;
  } else if (organ == "roots") {
    storage_capcity = 3.0;            // TODO: find a better value for this!
  } else if (organ == "all") {
    storage_capcity = 4.0;
  }

  if (nitrogen_balance > storage_capcity) {
    out = 1.0;
  } else {
    out = std::max(1.0 - pow(1.0 - nitrogen_balance/storage_capcity, 2.0), 0.0);
  }

  return(out);
}

double sugar_investment_for_nitrogen_target(double nitrogen_target,
                                            double ectomycorrhizal_uptake,
                                            double sugar_half_saturation) {
  // Ensure safe input values
  nitrogen_target = std::max(nitrogen_target, 0.0);
  ectomycorrhizal_uptake = std::max(ectomycorrhizal_uptake, 1e-12); // avoid divide-by-zero

  // Clamp nitrogen target just below max uptake
  double max_target = 0.999 * ectomycorrhizal_uptake;
  double safe_target = std::min(nitrogen_target, max_target);

  // Compute sugar investment
  return (sugar_half_saturation * safe_target) / (ectomycorrhizal_uptake - safe_target);
}

double nitrogen_transfer_from_sugar(double sugar_to_mycorrhiza,
                                    double ectomycorrhizal_uptake,
                                    double sugar_half_saturation) {
  // Clamp to avoid divide-by-zero
  sugar_to_mycorrhiza = std::max(sugar_to_mycorrhiza, 0.0);

  // Saturating function: more sugar → more N transferred
  double transfer_fraction = sugar_to_mycorrhiza / (sugar_to_mycorrhiza + sugar_half_saturation);

  // Total N delivered to tree
  return transfer_fraction * ectomycorrhizal_uptake;
}


uptake_structre nitrogen_uptake(double N,
                                double sugar_to_mycorrhiza,
                                double mycorrhizal_biomass,
                                double root_biomass,
                                double mycorrhizal_nitrogen_demand,
                                bool mycorrhiza_passive) {

  double uptake_capacity_constant = 0.1;
  double root_exploration = 0.0017;
  double mycelial_exploration = 0.5; // TODO: check numbers
  double collonialisation = 0.9; // TODO: don't make this dynamic yet! Just set values

  /*
   * Equations from Franklin 2014
   */

  double ectomycorrhizal_upatke_capacity = mycorrhizal_biomass * (uptake_capacity_constant * N * mycelial_exploration)/(uptake_capacity_constant + N * mycelial_exploration);
  double ectomycorrhizal_upatke = ectomycorrhizal_upatke_capacity * N / (ectomycorrhizal_upatke_capacity + N);

  double root_upatke_capacity = root_biomass * (uptake_capacity_constant * N * root_exploration)/(uptake_capacity_constant + N * root_exploration);
  double root_upatke = root_upatke_capacity * N / (root_upatke_capacity + N);

  double ectomycorrhizal_transfer;
  if (mycorrhiza_passive) {
    ectomycorrhizal_transfer = collonialisation * ectomycorrhizal_upatke;
  } else {
    ectomycorrhizal_transfer = nitrogen_transfer_from_sugar(sugar_to_mycorrhiza, ectomycorrhizal_upatke, 0.05); // TODO: tunable value
  }


  double total_uptake = root_upatke * (1 - collonialisation) + ectomycorrhizal_transfer;

  uptake_structre out;
  out.ectomycorrhizal_transfer = ectomycorrhizal_transfer;
  out.root_upatke = root_upatke;
  out.ectomycorrhizal_upatke = ectomycorrhizal_upatke;
  out.total_uptake = total_uptake;

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

double calculate_movement(double source_sugar, double source_starch, double source_min,
                          double source_capacity, double source_mass,
                          double target_sugar, double target_starch, double target_min,
                          double target_capacity, double target_mass) {
  if (source_mass <= 0.0 || target_mass <= 0.0) return 0.0;

  double source_relative = (source_sugar + source_starch - source_min) / (source_capacity * source_mass);
  double target_relative = (target_sugar + target_starch - target_min) / (target_capacity * target_mass);

  // Clamp both within [0,1]
  source_relative = std::min(std::max(source_relative, 0.0), 1.0);
  target_relative = std::min(std::max(target_relative, 0.0), 1.0);

  // Movement is the surplus in the source vs need in the target
  return std::max(source_relative - target_relative, 0.0);
}


/*
 * Some of the parameters are site dependent, could define the site side in the main function code so I just import the right variables into this code
 *
 * needles mass per year, should define the year before this function!
 */

carbo_balance sugar_model(int year,
                          int day,
                          double TAir,
                          double PAR,
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
                          double mycorrhizal_biomass,
                          carbo_tracker temperature_equilibrium, // Calculated in the main function

                          growth_out pot_growth,

                          carbo_tracker sugar,
                          carbo_tracker starch,

                          carbo_values_out parameters_in) {

  uptake_structre uptake;

  // TODO: phloem mass is hard to get within the model, maybe a product of diameter and height
  double phloem_mass = 7.410537931;
  double xylem_sh_mass = 74.10537931;
  double xylem_st_mass = 8.65862069;

  double phloem_respiration_share = phloem_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);
  double xylem_sh_respiration_share = xylem_sh_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);
  double xylem_st_respiration_share = xylem_st_mass / (phloem_mass + xylem_sh_mass + xylem_st_mass);

  double phloem_growth_share = phloem_mass / (phloem_mass + xylem_st_mass);
  double xylem_st_growth_share = xylem_st_mass / (phloem_mass + xylem_st_mass);

  // TODO: make sure that these make sense
  carbo_tracker storage_term;
  // 7.5% of mass should be for storage, von Arx 2017 max axis
  storage_term.needles  = storage_update_organs(resp.RmN + (1 + common.Rg_N) * (pot_growth.needles + pot_growth.bud), sugar.needles, parameters.lower_bound_needles, tree_alive);
  storage_term.phloem   = storage_update_organs(phloem_respiration_share * resp.RmS + phloem_growth_share * (1 + common.Rg_S) * (pot_growth.wall + pot_growth.height), sugar.phloem, parameters.lower_bound_phloem, tree_alive);
  storage_term.xylem_sh = storage_update_organs(xylem_sh_respiration_share * resp.RmS, sugar.xylem_sh, parameters.lower_bound_xylem_sh, tree_alive);
  storage_term.xylem_st = storage_update_organs(xylem_st_respiration_share * resp.RmS + xylem_st_growth_share * (1 + common.Rg_S) * (pot_growth.wall + pot_growth.height), sugar.xylem_st, parameters.lower_bound_xylem_st, tree_alive);
  storage_term.roots    = storage_update_organs(resp.RmR + (1 + common.Rg_R) * pot_growth.roots, sugar.roots, parameters.lower_bound_roots, tree_alive);

  growth_out nitrogen_capacity;
  nitrogen_capacity.needles = nitrogen_storage(nitrogen_balance/5, "needles");
  nitrogen_capacity.bud = nitrogen_storage(nitrogen_balance/5, "bud");
  nitrogen_capacity.wall = nitrogen_storage(nitrogen_balance/5, "wall");
  nitrogen_capacity.height = nitrogen_storage(nitrogen_balance/5, "height");
  nitrogen_capacity.roots = nitrogen_storage(nitrogen_balance/5, "roots");

  double growth = 0.0;
  double respiration_growth = 0.0;
  double respiration_maintainence = 0.0;

  // TOOD: clean the code!
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

      /*
       * Winter
       */

      carbo_tracker winter_costs;
      // 12% to 22% of water content reduces the freezing rate of water
      // Just using the values from autumn at the moment.
      // Temperature and PAR avaerage used to control the decrease / increase in allocation.

      // Lower bound added her so the processes don't happen unless a lower bound is met
      double winter = (1.0/(1.0 + std::exp(0.5*(TAir-5))) + 1.0/(1.0 + std::exp(0.2*(PAR-30))))/2.0;

      winter_costs.needles = 0.0536 * winter + parameters.lower_bound_needles;
      winter_costs.phloem = 0.2736 * winter + parameters.lower_bound_phloem;
      winter_costs.roots = 0.0300 * winter + parameters.lower_bound_roots;
      winter_costs.xylem_st = 0.0569 * winter + parameters.lower_bound_xylem_st;
      winter_costs.xylem_sh = 0.0231 * winter + parameters.lower_bound_xylem_sh;

      carbo_tracker storage_term;
      // 7.5% of mass should be for storage, von Arx 2017 max axis
      storage_term.needles  = storage_update_organs(resp.RmN + (1 + common.Rg_N) * (pot_growth.needles + pot_growth.bud), sugar.needles, winter_costs.needles, tree_alive);
      storage_term.phloem   = storage_update_organs(phloem_respiration_share * resp.RmS + phloem_growth_share * (1 + common.Rg_S) * (pot_growth.wall + pot_growth.height) + pot_growth.use - pot_growth.release, sugar.phloem, winter_costs.phloem, tree_alive);
      storage_term.xylem_sh = storage_update_organs(xylem_sh_respiration_share * resp.RmS, sugar.xylem_sh, winter_costs.xylem_sh, tree_alive);
      storage_term.xylem_st = storage_update_organs(xylem_st_respiration_share * resp.RmS + xylem_st_growth_share * (1 + common.Rg_S) * (pot_growth.wall + pot_growth.height), sugar.xylem_st, winter_costs.xylem_st, tree_alive);
      storage_term.roots    = storage_update_organs(resp.RmR + (1 + common.Rg_R) * pot_growth.roots, sugar.roots, winter_costs.roots, tree_alive);

      /*
       * Storage for respiration
       */

      carbo_tracker storage_term_resp;
      storage_term_resp.needles = respiration_check(storage_term.needles, 0.1); // TODO: randomly chose 0.1!
      storage_term_resp.phloem = respiration_check(storage_term.phloem, 0.1);
      storage_term_resp.xylem_sh = respiration_check(storage_term.xylem_sh, 0.1);
      storage_term_resp.xylem_st = respiration_check(storage_term.xylem_st, 0.1);
      storage_term_resp.roots = respiration_check(storage_term.roots, 0.1);

      /*
       * Balance calculations
       */

      double carbo_beginning = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + PF;
      double sugar_out_of_system = 0.0;

      /*
       * Growth and respiration
       */

      carbo_tracker growth_resp;
      growth_resp.needles =                               resp.RmN * storage_term_resp.needles  +                         (1 + common.Rg_N) * (std::min(storage_term.needles, nitrogen_capacity.needles)  * pot_growth.needles + std::min(storage_term.needles, nitrogen_capacity.bud) * pot_growth.bud);
      growth_resp.roots =                                 resp.RmR * storage_term_resp.roots    +                         (1 + common.Rg_R) * std::min(storage_term.roots, nitrogen_capacity.roots)       * pot_growth.roots;
      growth_resp.phloem =     phloem_respiration_share * resp.RmS * storage_term_resp.phloem   + phloem_growth_share   * (1 + common.Rg_S) * (std::min(storage_term.phloem, nitrogen_capacity.wall)      * pot_growth.wall    + std::min(storage_term.phloem, nitrogen_capacity.height) * pot_growth.height) + pot_growth.use - pot_growth.release;
      growth_resp.xylem_st = xylem_st_respiration_share * resp.RmS * storage_term_resp.xylem_st + xylem_st_growth_share * (1 + common.Rg_S) * (std::min(storage_term.xylem_st, nitrogen_capacity.wall)    * pot_growth.wall    + std::min(storage_term.xylem_st, nitrogen_capacity.height) * pot_growth.height);
      growth_resp.xylem_sh = xylem_sh_respiration_share * resp.RmS * storage_term_resp.xylem_sh;

      /*
       * Starch and sugar transfer
       */

      // Ad0_needles sugar to starch
      // delta_needles starch to sugar

      carbo_tracker to_starch;

      /*
       * Sugar and Starch Balence: Apply the transfer
       */

      if (sugar.needles < growth_resp.needles) {
        std::cout << "Day " << day + 1 << ": Stoage term is not capping growth enough needles\n";
      }
      sugar.needles = sugar.needles + PF - growth_resp.needles;


      if (sugar.phloem < growth_resp.phloem) {
        std::cout << "Day " << day + 1 << ": Stoage term is not capping growth enough phloem\n";
      }
      sugar.phloem = sugar.phloem - growth_resp.phloem;


      if (sugar.roots < growth_resp.roots) {
        std::cout << "Day " << day + 1 << ": Stoage term is not capping growth enough roots\n";
      }
      sugar.roots = sugar.roots - growth_resp.roots;


      if (sugar.xylem_sh < growth_resp.xylem_sh) {
        std::cout << "Day " << day + 1 << ": Stoage term is not capping growth enough xylem sh\n";
      }
      sugar.xylem_sh = sugar.xylem_sh - growth_resp.xylem_sh;


      if (sugar.xylem_st < growth_resp.xylem_st) {
        std::cout << "Day " << day + 1 << ": Stoage term is not capping growth enough xylem st\n";
      }
      sugar.xylem_st = sugar.xylem_st - growth_resp.xylem_st;

      /*
       * Starch
       */

      // === NEEDLES ===
      {
        double total = sugar.needles + starch.needles;
        double capacity = parameters.percentage_needle_storage * needles_mass;
        double space_left = std::max(capacity - total, 0.0);

        to_starch.needles = safe_to_starch_transfer(sugar.needles, winter_costs.needles, starch.needles,
                                                    parameters.Ad0_needles, parameters.delta_needles, space_left);

        sugar.needles -= to_starch.needles;
        starch.needles += to_starch.needles;

        if (sugar.needles < 0.0) std::cout << "Day " << day + 1 << " sugar.needles is negative after starch\n";
        if (starch.needles < 0.0) std::cout << "Day " << day + 1 << " starch.needles is negative after starch\n";
      }

      // === PHLOEM ===
      {
        double total = sugar.phloem + starch.phloem;
        double capacity = parameters.percentage_phloem_storage * phloem_mass;
        double space_left = std::max(capacity - total, 0.0);

        to_starch.phloem = safe_to_starch_transfer(sugar.phloem, winter_costs.phloem, starch.phloem,
                                                   parameters.Ad0_phloem, parameters.delta_phloem, space_left);

        sugar.phloem -= to_starch.phloem;
        starch.phloem += to_starch.phloem;

        if (sugar.phloem < 0.0) std::cout << "Day " << day + 1 << " sugar.phloem is negative after starch\n";
        if (starch.phloem < 0.0) std::cout << "Day " << day + 1 << " starch.phloem is negative after starch\n";
      }

      // === ROOTS ===
      {
        double total = sugar.roots + starch.roots;
        double capacity = parameters.percentage_roots_storage * root_mass;
        double space_left = std::max(capacity - total, 0.0);

        to_starch.roots = safe_to_starch_transfer(sugar.roots, winter_costs.roots, starch.roots,
                                                  parameters.Ad0_roots, parameters.delta_roots, space_left);

        sugar.roots -= to_starch.roots;
        starch.roots += to_starch.roots;

        if (sugar.roots < 0.0) std::cout << "Day " << day + 1 << " sugar.roots is negative after starch\n";
        if (starch.roots < 0.0) std::cout << "Day " << day + 1 << " starch.roots is negative after starch\n";
      }

      // === XYLEM SH ===
      {
        double total = sugar.xylem_sh + starch.xylem_sh;
        double capacity = parameters.percentage_xylem_sh_storage * xylem_sh_mass;
        double space_left = std::max(capacity - total, 0.0);

        to_starch.xylem_sh = safe_to_starch_transfer(sugar.xylem_sh, winter_costs.xylem_sh, starch.xylem_sh,
                                                     parameters.Ad0_xylem_sh, parameters.delta_xylem_sh, space_left);

        sugar.xylem_sh -= to_starch.xylem_sh;
        starch.xylem_sh += to_starch.xylem_sh;

        if (sugar.xylem_sh < 0.0) std::cout << "Day " << day + 1 << " sugar.xylem_sh is negative after starch\n";
        if (starch.xylem_sh < 0.0) std::cout << "Day " << day + 1 << " starch.xylem_sh is negative after starch\n";
      }

      // === XYLEM ST ===
      {
        double total = sugar.xylem_st + starch.xylem_st;
        double capacity = parameters.percentage_xylem_st_storage * xylem_st_mass;
        double space_left = std::max(capacity - total, 0.0);

        to_starch.xylem_st = safe_to_starch_transfer(sugar.xylem_st, winter_costs.xylem_st, starch.xylem_st,
                                                     parameters.Ad0_xylem_st, parameters.delta_xylem_st, space_left);

        sugar.xylem_st -= to_starch.xylem_st;
        starch.xylem_st += to_starch.xylem_st;

        if (sugar.xylem_st < 0.0) std::cout << "Day " << day + 1 << " sugar.xylem_st is negative after starch\n";
        if (starch.xylem_st < 0.0) std::cout << "Day " << day + 1 << " starch.xylem_st is negative after starch\n";
      }

      /*
       * Total actual growth
       */

      double growth_resp_combined = growth_resp.needles + growth_resp.phloem + growth_resp.xylem_sh + growth_resp.xylem_st + growth_resp.roots;

      sugar_out_of_system = growth_resp_combined;
      double carbo_ending = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + sugar_out_of_system;
      double difference = carbo_beginning - carbo_ending;
      if (difference > 1e-10) {
        std::cout << "On day " << day + 1 << " After starch transfer and growth: More carbohydrates compared to the BEGINNING of the day. Difference is: " << difference << "\n";
      } else if (difference < -1e-10) {
        std::cout << "On day " << day + 1 << " After starch transfer and growth: More carbohydrate compared to the END of the day. Difference is: " << difference << "\n";
      }

      /*
       * The concentration gradient is represented by the difference in storage terms.
       *
       * This is then multiplied by the total storage in the organ and divided by a delay factor.
       *
       * The amount  that is transferred should then either be the sugar that is over the storage capacity of the organ
       */

      /*
       * === Mycorrhizal sugar export ===
       */

      // --- 1. Mycorrhizal demand based on nutrient contrast ---
      // Step 1: Compute current total sugar + starch
      double total_sugar_now = sugar.needles + sugar.phloem + sugar.roots + sugar.xylem_sh + sugar.xylem_st;
      double total_starch_now = starch.needles + starch.phloem + starch.roots + starch.xylem_sh + starch.xylem_st;
      double total_storage = total_sugar_now + total_starch_now;

      double myco_demand = parameters.mycorrhiza_threshold * total_sugar_now;

      if (surplus_c) {
        myco_demand = 0.0;
      } else if (nitrogen_contrast) {
        // Step 1: Compute nitrogen deficit (targeted increase)
        double nitrogen_capacity_average = (
          nitrogen_capacity.needles +
            nitrogen_capacity.bud +
            nitrogen_capacity.wall +
            nitrogen_capacity.height +
            nitrogen_capacity.roots
        ) / 5.0;

        double average_storage_term = (
          storage_term.needles +
            storage_term.phloem +
            storage_term.roots +
            storage_term.xylem_sh +
            storage_term.xylem_st
        ) / 5.0;

        double normalised_nitrogen_target_increase = std::max(average_storage_term - nitrogen_capacity_average, 0.0);

        // Step 2: Estimate mycorrhizal uptake capacity at current environment
        double N = 0.5;
        double ecto_uptake_capacity = mycorrhizal_biomass * (
          (0.1 * N * 0.5) / (0.1 + N * 0.5)  // from Franklin 2014, with uptake_capacity_constant = 0.1, mycelial_exploration = 0.5
        );
        double ectomycorrhizal_uptake = ecto_uptake_capacity * N / (ecto_uptake_capacity + N);

        // Step 3: Estimate sugar investment required to obtain that N increase
        double sugar_half_saturation = 0.05;  // tunable parameter
        double nitrogen_target = normalised_nitrogen_target_increase;

        double sugar_required = sugar_investment_for_nitrogen_target(
          nitrogen_target,
          ectomycorrhizal_uptake,
          sugar_half_saturation
        );

        // Step 4: Cap investment by available sugar budget
        double total_use = growth_resp.needles + growth_resp.phloem + growth_resp.xylem_sh + growth_resp.xylem_st + growth_resp.roots;

        double sugar_extra = std::max(total_storage - total_use, 0.0);

        // Step 5: Final sugar allocated to mycorrhiza, capped by what's available
        myco_demand = std::min(sugar_required, sugar_extra);
      }


      /*
       * === Excess sugar export beyond capacity ===
       */

      // check the mins and max as well as the order here
      // mini functions?

      total_sugar_now -= myco_demand;

      // Step 2: Compute total organ storage capacity
      double total_capacity =
        parameters.percentage_needle_storage * needles_mass +
        parameters.percentage_phloem_storage * phloem_mass +
        parameters.percentage_roots_storage * root_mass +
        parameters.percentage_xylem_sh_storage * xylem_sh_mass +
        parameters.percentage_xylem_st_storage * xylem_st_mass;

      // Step 3: Calculate overflow sugar to export
      double excess_carbo = std::max(total_storage - total_capacity, 0.0);

      // Step 4: Apply export rate
      double export_rate = 0.5;  // Reuse this as export efficiency
      double export_amount = export_rate * excess_carbo;

      // Step 5: Transfer excess sugar into roots/phloem first (push through tree)

      // First from phloem → roots
      double root_capacity = parameters.percentage_roots_storage * root_mass;
      double root_current = sugar.roots + starch.roots;
      double root_space = std::max(root_capacity - root_current, 0.0);

      double transfer_to_roots = std::min({sugar.phloem, root_space, export_amount});
      sugar.phloem -= transfer_to_roots;
      sugar.roots += transfer_to_roots;
      export_amount -= transfer_to_roots;

      // Step 6: Export from roots (main interface)
      double export_from_roots = std::min(sugar.roots, export_amount);
      sugar.roots -= export_from_roots;
      export_amount -= export_from_roots;

      // Step 7: Export remainder from phloem (if roots not enough)
      double export_from_phloem = std::min(sugar.phloem, export_amount);
      sugar.phloem -= export_from_phloem;
      export_amount -= export_from_phloem;

      // Final export sum
      double total_exported = export_from_roots + export_from_phloem;
      sugar_out_of_system += total_exported;

      /*
       * Concentration
       */

      // --- INPUTS ---
      double sugar_needles = sugar.needles;
      double sugar_phloem  = sugar.phloem;

      // --- STEP 1: Compute downstream demand ---
      double demand_roots    = growth_resp.roots + winter_costs.needles + myco_demand;
      double demand_xylem_st = growth_resp.xylem_st + winter_costs.xylem_st;
      double demand_xylem_sh = growth_resp.xylem_sh + winter_costs.xylem_sh;

      double total_demand = demand_roots + demand_xylem_st + demand_xylem_sh;

      // --- STEP 2: Determine how much sugar phloem still needs ---
      double available_phloem = sugar_phloem - winter_costs.phloem;
      double phloem_to_cover_sinks = std::max(total_demand - available_phloem, 0.0);

      // --- STEP 3: Cap needle export by available sugar ---
      double needle_reserve = std::max(sugar_needles - winter_costs.needles, 0.0);
      double needle_to_phloem = std::min(needle_reserve, phloem_to_cover_sinks);

      conc_gradient concentration_gradient;
      concentration_gradient.needles_to_phloem = needle_to_phloem;

      // Apply transfer
      sugar.needles -= needle_to_phloem;
      sugar.phloem  += needle_to_phloem;

      if (sugar.needles < 0.0) std::cout << "sugar.needles is negative after needle concentration!\n";
      if (sugar.phloem < 0.0) std::cout << "sugar.needles is negative after needle concentration!\n";

      // STEP 4: Redistribute sugar in phloem to sinks, capped by supply but keep minimum reserve
      // STEP 4: Redistribute sugar in phloem to sinks, constrained by demand, capacity, and supply

      double phloem_min_reserve = winter_costs.phloem;
      double phloem_available_for_transfer = std::max(sugar.phloem - phloem_min_reserve, 0.0);

      // === ROOTS ===
      double root_storage_total = sugar.roots + starch.roots;
      double root_capacity_max = parameters.percentage_roots_storage * root_mass;
      double root_capacity_remaining = std::max(root_capacity_max - root_storage_total, 0.0);

      // Only growth + winter demand here; myco_demand is handled separately
      double root_growth_and_maintenance_demand = growth_resp.roots + winter_costs.roots;

      double root_transfer = std::min({phloem_available_for_transfer, root_growth_and_maintenance_demand, root_capacity_remaining});
      sugar.phloem -= root_transfer;
      sugar.roots  += root_transfer;
      phloem_available_for_transfer -= root_transfer;
      concentration_gradient.phloem_to_roots = root_transfer;

      // == MYCORRHIZAL EXPORT ==

      // Only export based on what's actually in roots now, capped by demand
      double myco_export = std::min(sugar.roots, myco_demand);
      sugar.roots -= myco_export;
      concentration_gradient.roots_to_myco = myco_export;
      sugar_out_of_system += myco_export;

      // === XYLEM ST ===
      double xylem_st_storage_total = sugar.xylem_st + starch.xylem_st;
      double xylem_st_capacity_max = parameters.percentage_xylem_st_storage * xylem_st_mass;
      double xylem_st_capacity_remaining = std::max(xylem_st_capacity_max - xylem_st_storage_total, 0.0);

      double xylem_st_demand = demand_xylem_st;
      double xylem_st_transfer = std::min({phloem_available_for_transfer, xylem_st_demand, xylem_st_capacity_remaining});
      sugar.phloem -= xylem_st_transfer;
      sugar.xylem_st += xylem_st_transfer;
      phloem_available_for_transfer -= xylem_st_transfer;
      concentration_gradient.phloem_to_xylem_st = xylem_st_transfer;

      // === XYLEM SH ===
      double xylem_sh_storage_total = sugar.xylem_sh + starch.xylem_sh;
      double xylem_sh_capacity_max = parameters.percentage_xylem_sh_storage * xylem_sh_mass;
      double xylem_sh_capacity_remaining = std::max(xylem_sh_capacity_max - xylem_sh_storage_total, 0.0);

      double xylem_sh_demand = demand_xylem_sh;
      double xylem_sh_transfer = std::min({phloem_available_for_transfer, xylem_sh_demand, xylem_sh_capacity_remaining});
      sugar.phloem -= xylem_sh_transfer;
      sugar.xylem_sh += xylem_sh_transfer;
      phloem_available_for_transfer -= xylem_sh_transfer;
      concentration_gradient.phloem_to_xylem_sh = xylem_sh_transfer;

      // === SAFETY CHECKS ===
      if (sugar.phloem < 0.0) std::cout << "sugar.phloem is negative after all concentration!\n";
      if (sugar.roots < 0.0) std::cout << "sugar.roots is negative after all concentration!\n";
      if (sugar.xylem_st < 0.0) std::cout << "sugar.xylem_st is negative after all concentration!\n";
      if (sugar.xylem_sh < 0.0) std::cout << "sugar.xylem_sh is negative after all concentration!\n";

      /*
       * Concentration
       */

      double carbo_ending_2 = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + sugar_out_of_system;
      double difference_2 = carbo_beginning - carbo_ending_2;
      if (difference_2 > 1e-8) {
        std::cout << "On day " << day + 1 << " After concentration: More carbohydrates compared to the BEGINNING of the day. Difference is: " << difference_2 << "\n";
      } else if (difference_2 < -1e-8) {
        std::cout << "On day " << day + 1 << " After concentration: More carbohydrate compared to the END of the day. Difference is: " << difference_2 << "\n";
      }

      /*
       * Nitrogen
       */

      if (nitrogen_change) {
        nitrogen_capacity.needles = nitrogen_storage(nitrogen_balance/5.0, "needles");
        nitrogen_capacity.bud     = nitrogen_storage(nitrogen_balance/5.0, "bud");
        nitrogen_capacity.wall    = nitrogen_storage(nitrogen_balance/5.0, "wall");
        nitrogen_capacity.height  = nitrogen_storage(nitrogen_balance/5.0, "height");
        nitrogen_capacity.roots   = nitrogen_storage(nitrogen_balance/5.0, "roots");

        double N = 0.5;

        double mycorrhizal_nitrogen_demand = 0.2; // NOTE: currently not used
        uptake = nitrogen_uptake(N,
                                 myco_demand,
                                 mycorrhizal_biomass,
                                 root_mass,
                                 mycorrhizal_nitrogen_demand,
                                 FALSE);

        // C:N ratios are from the Korhonen 2013 paper
        // TODO: real value for the roots
        // TODO: I think the phloem share should be included here!
        nitrogen_balance = nitrogen_balance + uptake.total_uptake - growth_resp.needles - growth_resp.phloem - growth_resp.xylem_sh - growth_resp.xylem_st - growth_resp.roots;
      }

      /*
       * Storage check
       */

      if ((starch.needles == 0) && (starch.phloem == 0) && (starch.xylem_sh == 0) && (starch.xylem_st == 0) && (starch.roots == 0) &&
          (sugar.needles == 0) && (sugar.phloem == 0) && (sugar.xylem_sh == 0) && (sugar.xylem_st == 0) && (sugar.roots == 0)) {
        std::cerr << " Day " << day << " No total storage - plant died!" << "\n";
        // tree_alive = false;
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
    double ak = 1 / (1 - 1/exp(parameters.alfa * (parameters.sugar00 + parameters.starch00 - parameters.Wala)));
    double storage, storage_term_Rm, sugar_all, starch_all, to_sugar, to_starch;
    double nitrogen_capacity_all;
    double myco_allocation;
    if (day == 0) {
      sugar_all = parameters.sugar0;
      starch_all = parameters.starch0;
      to_sugar = 0;
      to_starch = 0;
      storage = storage_term.respiration = 1;
      nitrogen_capacity_all = 1;
    } else {
      storage = std::max(0.0 , std::min(1.0 , ak * (1.0 - 1.0 / exp(parameters.alfa * (sugar.needles + starch.needles - parameters.Wala)))));
      nitrogen_capacity_all = nitrogen_storage(nitrogen_balance/5.0, "all");
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

      sugar_all = sugar.needles + PF - pot_growth.use + pot_growth.release - std::min(storage_term.respiration, nitrogen_capacity_all) * resp.Rm_a -
        (1 + common.Rg_S) * std::min(storage, nitrogen_capacity_all) * pot_growth.height -
        (1 + common.Rg_S) * std::min(storage, nitrogen_capacity_all) * pot_growth.diameter -
        (1 + common.Rg_N) * std::min(storage, nitrogen_capacity_all) * pot_growth.needles -
        (1 + common.Rg_R) * std::min(storage, nitrogen_capacity_all) * pot_growth.roots -
        (1 + common.Rg_N) * std::min(storage, nitrogen_capacity_all) * pot_growth.bud -
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

    nitrogen_capacity.needles = nitrogen_capacity_all;
    nitrogen_capacity.bud = nitrogen_capacity_all;
    nitrogen_capacity.roots = nitrogen_capacity_all;
    nitrogen_capacity.wall = nitrogen_capacity_all;
    nitrogen_capacity.height = nitrogen_capacity_all;

    respiration_growth = (common.Rg_S) * std::min(storage, nitrogen_capacity_all) * pot_growth.height -
      (common.Rg_S) * std::min(storage, nitrogen_capacity_all) * pot_growth.diameter -
      (common.Rg_N) * std::min(storage, nitrogen_capacity_all) * pot_growth.needles -
      (common.Rg_R) * std::min(storage, nitrogen_capacity_all) * pot_growth.roots -
      (common.Rg_N) * std::min(storage, nitrogen_capacity_all) * pot_growth.bud;
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

  if (nitrogen_change) {
    out.uptake = uptake;
  }

  return out;
}
