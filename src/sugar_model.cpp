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


double sigma(double input, double zeta) {
  double out = 1.0/(1.0 + std::exp(-zeta*input));
  return out;
}

double safe_to_starch_transfer(
    double sugar,
    double winter,
    double starch,
    double Ad0,
    double delta,
    double max_capacity
) {
  double zeta = 20;

  // Sugar surplus after winter needs
  double S = sugar - winter;

  double sugar_limitation = sigma(S, zeta);
  double starch_limitation = sigma(starch, zeta);

  double exponent = max_capacity - S - starch;
  double max_capacity_limitation = sigma(exponent, zeta);

  double synthesis = Ad0 * S * max_capacity_limitation;
  double remobilisation = delta * starch * starch_limitation;

  double proposed_flux = sugar_limitation * synthesis - (1 - sugar_limitation) * remobilisation;

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
    out = std::min(nitrogen_balance/storage_capcity, 1.0);
  }

  return(out);
}

double nitrogen_transfer_from_mycorrhiza(double sugar_to_mycorrhiza,
                                    double sugar_half_saturation) {
  // Clamp to avoid divide-by-zero
  if (sugar_to_mycorrhiza < 0.0) {
    std::cout << "Warning: sugar_to_mycorrhiza is less than zero, claped to avoid erros.";
  }
  sugar_to_mycorrhiza = std::max(sugar_to_mycorrhiza, 0.0);

  // Saturating function: more sugar → more N transferred
  double transfer_fraction = sugar_to_mycorrhiza / (sugar_to_mycorrhiza + sugar_half_saturation);

  // Total N delivered to tree
  return transfer_fraction;
}


uptake_structre nitrogen_uptake(double N,
                                double sugar_to_mycorrhiza,
                                double surplus_sugar,
                                double mycorrhizal_biomass,
                                double root_biomass,
                                double mycorrhizal_nitrogen_demand,
                                bool mycorrhiza_passive) {

  double uptake_capacity_constant = 0.1;
  double root_exploration = 0.0017;
  double mycelial_exploration = 0.5; // TODO: check numbers

  /*
   * Equations from Franklin 2014
   */

  double ectomycorrhizal_uptake_capacity = mycorrhizal_biomass * (uptake_capacity_constant * N * mycelial_exploration)/(uptake_capacity_constant + N * mycelial_exploration);
  double ectomycorrhizal_uptake = ectomycorrhizal_uptake_capacity * N / (ectomycorrhizal_uptake_capacity + N);

  double root_uptake_capacity = root_biomass * (uptake_capacity_constant * N * root_exploration)/(uptake_capacity_constant + N * root_exploration);
  double root_uptake = root_uptake_capacity * N / (root_uptake_capacity + N);

  double ectomycorrhizal_transfer;
  if (mycorrhiza_passive) {
    ectomycorrhizal_transfer = nitrogen_transfer_from_mycorrhiza(sugar_to_mycorrhiza, 0.00005) * ectomycorrhizal_uptake;
  } else {
    std::cout << "Write this :)\n";
  }

  double free_sugar_percentage = surplus_sugar/(surplus_sugar + sugar_to_mycorrhiza + 1e-10);

  double total_uptake = ectomycorrhizal_transfer + free_sugar_percentage * root_uptake;

  uptake_structre out;
  out.ectomycorrhizal_transfer = ectomycorrhizal_transfer;
  out.root_uptake = root_uptake;
  out.ectomycorrhizal_uptake = ectomycorrhizal_uptake;
  out.total_uptake = total_uptake;

  return(out);
}

double costant_excess(double nitrogen_target, double sugar_to_invest) {

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

void sugar_model(int year,
                 int days_gone,
                 int day,
                 double TAir,
                 double PAR,
                 double PF,

                 const CASSIA_common& common,
                 const CASSIA_parameters& parameters,

                 double D00,
                 const growth_state& tree_state,

                 double& nitrogen_balance,
                 bool nitrogen_change,
                 bool nitrogen_contrast,

                 Settings boolsettings,

                 bool& tree_alive,
                 bool surplus_c,

                 double needles_mass, // Repola

                 carbo_tracker& sugar,
                 carbo_tracker& starch,
                 carbo_tracker& storage_term,
                 growth_out& nitrogen_capacity,
                 output_vector& out) {

  uptake_structre uptake;

  // TODO: phloem mass is hard to get within the model, maybe a product of diameter and height
  // double phloem_mass = 7.410537931;
  // double xylem_sh_mass = 74.10537931;
  double xylem_st_mass = 8.65862069;

  double phloem_respiration_share = out.culm_growth.phloem[day-1] / (out.culm_growth.phloem[day-1] + out.culm_growth.xylem_sh[day-1] + xylem_st_mass);
  double xylem_sh_respiration_share = out.culm_growth.xylem_sh[day-1] / (out.culm_growth.phloem[day-1] + out.culm_growth.xylem_sh[day-1] + xylem_st_mass);
  double xylem_st_respiration_share = xylem_st_mass / (out.culm_growth.phloem[day-1] + out.culm_growth.xylem_sh[day-1] + xylem_st_mass);

  double phloem_growth_share = out.culm_growth.phloem[day-1] / (out.culm_growth.phloem[day-1] + out.culm_growth.xylem_sh[day-1]);
  double xylem_sh_growth_share = out.culm_growth.xylem_sh[day-1] / (out.culm_growth.phloem[day-1] + out.culm_growth.xylem_sh[day-1]);

  // 7.5% of mass should be for storage, von Arx 2017 max axis
  storage_term.needles  = storage_update_organs(tree_state.RmN + (1 + common.Rg_N) * (tree_state.needles + tree_state.bud), sugar.needles, parameters.lower_bound_needles, tree_alive);
  storage_term.phloem   = storage_update_organs(phloem_respiration_share * tree_state.RmS + phloem_growth_share * (1 + common.Rg_S) * (tree_state.wall + tree_state.height), sugar.phloem, parameters.lower_bound_phloem, tree_alive);
  storage_term.xylem_sh = storage_update_organs(xylem_sh_respiration_share * tree_state.RmS + xylem_sh_growth_share * (1 + common.Rg_S) * (tree_state.wall + tree_state.height), sugar.xylem_sh, parameters.lower_bound_xylem_sh, tree_alive);
  storage_term.xylem_st = storage_update_organs(xylem_st_respiration_share * tree_state.RmS, sugar.xylem_st, parameters.lower_bound_xylem_st, tree_alive);
  storage_term.roots    = storage_update_organs(tree_state.RmR + (1 + common.Rg_R) * tree_state.roots, sugar.roots, parameters.lower_bound_roots, tree_alive);

  nitrogen_capacity.needles = nitrogen_storage(nitrogen_balance/5, "needles");
  nitrogen_capacity.bud = nitrogen_storage(nitrogen_balance/5, "bud");
  nitrogen_capacity.wall = nitrogen_storage(nitrogen_balance/5, "wall");
  nitrogen_capacity.height = nitrogen_storage(nitrogen_balance/5, "height");
  nitrogen_capacity.roots = nitrogen_storage(nitrogen_balance/5, "roots");

  double growth = 0.0;
  // TODO: I don't think I calculate the below anymore
  double respiration_growth = 0.0;
  double respiration_maintainence = 0.0;

  // TOOD: clean the code!
  double sB0;
  carbo_tracker As_new, Ad_new;

  if (boolsettings.sperling_model) {
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

      // 7.5% of mass should be for storage, von Arx 2017 max axis
      storage_term.needles  = storage_update_organs(tree_state.RmN + (1 + common.Rg_N) * (tree_state.needles + tree_state.bud), sugar.needles, winter_costs.needles, tree_alive);
      storage_term.phloem   = storage_update_organs(phloem_respiration_share * tree_state.RmS + phloem_growth_share * (1 + common.Rg_S) * (tree_state.wall + tree_state.height) + tree_state.use - tree_state.release, sugar.phloem, winter_costs.phloem, tree_alive);
      storage_term.xylem_sh = storage_update_organs(xylem_sh_respiration_share * tree_state.RmS + xylem_sh_growth_share * (1 + common.Rg_S) * (tree_state.wall + tree_state.height), sugar.xylem_sh, winter_costs.xylem_sh, tree_alive);
      storage_term.xylem_st = storage_update_organs(xylem_st_respiration_share * tree_state.RmS, sugar.xylem_st, winter_costs.xylem_st, tree_alive);
      storage_term.roots    = storage_update_organs(tree_state.RmR + (1 + common.Rg_R) * tree_state.roots, sugar.roots, winter_costs.roots, tree_alive);

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
      growth_resp.needles =                               tree_state.RmN * storage_term_resp.needles  +                         (1 + common.Rg_N) * (std::min(storage_term.needles, nitrogen_capacity.needles)  * tree_state.needles + std::min(storage_term.needles, nitrogen_capacity.bud) * tree_state.bud);
      growth_resp.roots =                                 tree_state.RmR * storage_term_resp.roots    +                         (1 + common.Rg_R) * std::min(storage_term.roots, nitrogen_capacity.roots)       * tree_state.roots;
      growth_resp.phloem =     phloem_respiration_share * tree_state.RmS * storage_term_resp.phloem   + phloem_growth_share   * (1 + common.Rg_S) * (std::min(storage_term.phloem, nitrogen_capacity.wall)      * tree_state.wall    + std::min(storage_term.phloem, nitrogen_capacity.height) * tree_state.height) + tree_state.use - tree_state.release;
      growth_resp.xylem_st = xylem_st_respiration_share * tree_state.RmS * storage_term_resp.xylem_st;
      growth_resp.xylem_sh = xylem_sh_respiration_share * tree_state.RmS * storage_term_resp.xylem_sh + xylem_sh_growth_share * (1 + common.Rg_S) * (std::min(storage_term.xylem_st, nitrogen_capacity.wall)    * tree_state.wall    + std::min(storage_term.xylem_st, nitrogen_capacity.height) * tree_state.height);

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
        // Calculate the difference
        double growth_difficit = storage_term.needles - nitrogen_capacity.needles +
          storage_term.needles - nitrogen_capacity.bud +
          storage_term.roots - nitrogen_capacity.roots +
          phloem_growth_share * storage_term.phloem + xylem_sh_growth_share * storage_term.xylem_sh - nitrogen_capacity.wall +
          phloem_growth_share * storage_term.phloem + xylem_sh_growth_share * storage_term.xylem_sh - nitrogen_capacity.height;

        if (growth_difficit <= 0.0) {
          // More sugar, so sugar is limiting the growth
          myco_demand = 0.0;
        } else if (growth_difficit > 0.0) {
          // TODO: share of the organs correct
          double nitrogen_target = 4.9/0.5 * (storage_term.needles - nitrogen_capacity.needles) * tree_state.needles -
            2.1/0.5 * (storage_term.needles - nitrogen_capacity.bud) * tree_state.bud -
            3.7/0.5 * phloem_growth_share * ((storage_term.phloem - nitrogen_capacity.wall) * tree_state.wall -
            4.9/0.5 * (storage_term.phloem - nitrogen_capacity.height) * tree_state.height) -
            3.7/0.5 * xylem_sh_growth_share * ((storage_term.xylem_st - nitrogen_capacity.wall) * tree_state.wall -
            4.9/0.5 * (storage_term.xylem_st - nitrogen_capacity.height) * tree_state.height) -
            3.0/0.5 * (storage_term.roots - nitrogen_capacity.roots) * tree_state.roots;

          double sugar_to_invest = (1 + common.Rg_N) * ((storage_term.needles - nitrogen_capacity.needles) * tree_state.bud + (storage_term.needles - nitrogen_capacity.bud) * tree_state.bud) +
            (1 + common.Rg_R) * (storage_term.roots - nitrogen_capacity.roots) * tree_state.roots +
            phloem_growth_share * (1 + common.Rg_S) * ((storage_term.phloem - nitrogen_capacity.wall) * tree_state.wall + (storage_term.phloem - nitrogen_capacity.height) * tree_state.height) +
            xylem_sh_growth_share * (1 + common.Rg_S) * ((storage_term.xylem_st - nitrogen_capacity.wall) * tree_state.wall + (storage_term.xylem_st - nitrogen_capacity.height) * tree_state.height);

          myco_demand = costant_excess(nitrogen_target, sugar_to_invest);
        }
      }

      /*
       * Concentration and sugar redistribution
       */

      // --- INPUTS ---
      double sugar_needles = sugar.needles;
      double sugar_phloem  = sugar.phloem;

      // --- STEP 1: Compute sink demands and capacities ---

      // ROOTS
      double root_storage_total = sugar.roots + starch.roots;
      double root_capacity_max = parameters.percentage_roots_storage * out.culm_growth.roots[day-1];
      double root_capacity_remaining = std::max(root_capacity_max - root_storage_total, 0.0);
      // TODO: is the winter cost double counted here?
      double demand_roots = growth_resp.roots + winter_costs.roots;
      double actual_root_demand = std::min(demand_roots, root_capacity_remaining);

      // XYLEM ST
      double xylem_st_storage_total = sugar.xylem_st + starch.xylem_st;
      double xylem_st_capacity_max = parameters.percentage_xylem_st_storage * xylem_st_mass;
      double xylem_st_capacity_remaining = std::max(xylem_st_capacity_max - xylem_st_storage_total, 0.0);
      double demand_xylem_st = growth_resp.xylem_st + winter_costs.xylem_st;
      double actual_xylem_st_demand = std::min(demand_xylem_st, xylem_st_capacity_remaining);

      // XYLEM SH
      double xylem_sh_storage_total = sugar.xylem_sh + starch.xylem_sh;
      double xylem_sh_capacity_max = parameters.percentage_xylem_sh_storage * out.culm_growth.xylem_sh[day-1];
      double xylem_sh_capacity_remaining = std::max(xylem_sh_capacity_max - xylem_sh_storage_total, 0.0);
      double demand_xylem_sh = growth_resp.xylem_sh + winter_costs.xylem_sh;
      double actual_xylem_sh_demand = std::min(demand_xylem_sh, xylem_sh_capacity_remaining);

      // PHLOEM SELF-DEMAND
      double demand_phloem = growth_resp.phloem + winter_costs.phloem;

      // Compute total effective (capacity-limited + flow-through + phloem) demand
      double capacity_limited_demand = actual_root_demand + actual_xylem_st_demand + actual_xylem_sh_demand;
      double total_effective_demand = capacity_limited_demand + myco_demand + demand_phloem;

      // --- STEP 2: Compute phloem gap to meet demands ---
      double available_phloem = sugar_phloem;
      double phloem_to_cover_sinks = std::max(total_effective_demand - available_phloem, 0.0);

      // --- STEP 3: Determine needle-to-phloem transfer (capped by reserve and phloem need) ---
      // TODO: is there a double counting of the growth_resp here?Need
      double needle_reserve = std::max(sugar_needles - winter_costs.needles - growth_resp.needles, 0.0);
      double needle_to_phloem = std::min(needle_reserve, phloem_to_cover_sinks);

      // Apply transfer
      conc_gradient concentration_gradient;
      concentration_gradient.needles_to_phloem = needle_to_phloem;

      sugar.needles -= needle_to_phloem;
      sugar.phloem  += needle_to_phloem;

      if (sugar.needles < 0.0) std::cout << "sugar.needles is negative after needle concentration!\n";
      if (sugar.phloem  < 0.0) std::cout << "sugar.phloem is negative after needle concentration!\n";

      // === STEP 4: Redistribute sugar in phloem to sinks (roots, mycorrhiza, xylem) ===

      // Maintain reserve in phloem for its own costs
      double phloem_min_reserve = demand_phloem;
      double phloem_available_for_transfer = std::max(sugar.phloem - phloem_min_reserve, 0.0);

      // === ROOTS ===
      double root_growth_and_maintenance_demand = demand_roots;
      // === ROOTS ===
      double root_storage_transfer = std::min(phloem_available_for_transfer, actual_root_demand);
      sugar.phloem -= root_storage_transfer;
      sugar.roots += root_storage_transfer;
      phloem_available_for_transfer -= root_storage_transfer;
      concentration_gradient.phloem_to_roots = root_storage_transfer;

      // === XYLEM ST ===
      double xylem_st_transfer = std::min({phloem_available_for_transfer, actual_xylem_st_demand});
      sugar.phloem -= xylem_st_transfer;
      sugar.xylem_st += xylem_st_transfer;
      phloem_available_for_transfer -= xylem_st_transfer;
      concentration_gradient.phloem_to_xylem_st = xylem_st_transfer;

      // === XYLEM SH ===
      double xylem_sh_transfer = std::min({phloem_available_for_transfer, actual_xylem_sh_demand});
      sugar.phloem -= xylem_sh_transfer;
      sugar.xylem_sh += xylem_sh_transfer;
      phloem_available_for_transfer -= xylem_sh_transfer;
      concentration_gradient.phloem_to_xylem_sh = xylem_sh_transfer;

      // === MYCORRHIZA: Direct export from phloem ===
      double myco_export = std::min(phloem_available_for_transfer, myco_demand);
      sugar.phloem -= myco_export;
      phloem_available_for_transfer -= myco_export;

      concentration_gradient.roots_to_myco = myco_export;
      sugar_out_of_system += myco_export;
      sugar.mycorrhiza = myco_export;

      // --- BALANCE CHECK ---
      double carbo_ending_2 = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots +
        sugar_out_of_system;

      double difference_2 = carbo_beginning - carbo_ending_2;

      if (difference_2 > 1e-8) {
        std::cout << "Day " << day + 1 << ": More carbs at beginning. Difference: " << difference_2 << "\n";
      } else if (difference_2 < -1e-8) {
        std::cout << "Day " << day + 1 << ": More carbs at end. Difference: " << difference_2 << "\n";
      }

      /*
       * === Starch and Sugar ===
       */

      // === NEEDLES ===
      {
        double capacity = parameters.percentage_needle_storage * needles_mass;

        to_starch.needles = safe_to_starch_transfer(
          sugar.needles, winter_costs.needles, starch.needles,
          parameters.Ad0_needles, parameters.delta_needles, capacity
        );

        sugar.needles -= to_starch.needles;
        starch.needles += to_starch.needles;

        if (sugar.needles < 0.0) std::cout << "Day " << day + 1 << " sugar.needles is negative after starch\n";
        if (starch.needles < 0.0) std::cout << "Day " << day + 1 << " starch.needles is negative after starch\n";
      }

      // === PHLOEM ===
      {
        double capacity = parameters.percentage_phloem_storage * out.culm_growth.phloem[day-1];

        to_starch.phloem = safe_to_starch_transfer(
          sugar.phloem, winter_costs.phloem, starch.phloem,
          parameters.Ad0_phloem, parameters.delta_phloem, capacity
        );

        sugar.phloem -= to_starch.phloem;
        starch.phloem += to_starch.phloem;

        if (sugar.phloem < 0.0) std::cout << "Day " << day + 1 << " sugar.phloem is negative after starch\n";
        if (starch.phloem < 0.0) std::cout << "Day " << day + 1 << " starch.phloem is negative after starch\n";
      }

      // === ROOTS ===
      {
        double capacity = parameters.percentage_roots_storage * out.culm_growth.roots[day-1];

        to_starch.roots = safe_to_starch_transfer(
          sugar.roots, winter_costs.roots, starch.roots,
          parameters.Ad0_roots, parameters.delta_roots, capacity
        );

        sugar.roots -= to_starch.roots;
        starch.roots += to_starch.roots;

        if (sugar.roots < 0.0) std::cout << "Day " << day + 1 << " sugar.roots is negative after starch\n";
        if (starch.roots < 0.0) std::cout << "Day " << day + 1 << " starch.roots is negative after starch\n";
      }

      // === XYLEM SH ===
      {
        double capacity = parameters.percentage_xylem_sh_storage * out.culm_growth.xylem_sh[day-1];

        to_starch.xylem_sh = safe_to_starch_transfer(
          sugar.xylem_sh, winter_costs.xylem_sh, starch.xylem_sh,
          parameters.Ad0_xylem_sh, parameters.delta_xylem_sh, capacity
        );

        sugar.xylem_sh -= to_starch.xylem_sh;
        starch.xylem_sh += to_starch.xylem_sh;

        if (sugar.xylem_sh < 0.0) std::cout << "Day " << day + 1 << " sugar.xylem_sh is negative after starch\n";
        if (starch.xylem_sh < 0.0) std::cout << "Day " << day + 1 << " starch.xylem_sh is negative after starch\n";
      }

      // === XYLEM ST ===
      {
        double capacity = parameters.percentage_xylem_st_storage * xylem_st_mass;

        to_starch.xylem_st = safe_to_starch_transfer(
          sugar.xylem_st, winter_costs.xylem_st, starch.xylem_st,
          parameters.Ad0_xylem_st, parameters.delta_xylem_st, capacity
        );

        sugar.xylem_st -= to_starch.xylem_st;
        starch.xylem_st += to_starch.xylem_st;

        if (sugar.xylem_st < 0.0) std::cout << "Day " << day + 1 << " sugar.xylem_st is negative after starch\n";
        if (starch.xylem_st < 0.0) std::cout << "Day " << day + 1 << " starch.xylem_st is negative after starch\n";
      }

      /*
       * === Excess sugar export beyond capacity ===
       */

      // --- STEP 1: Compute per-organ storage capacities ---
      double needle_capacity = parameters.percentage_needle_storage * needles_mass;
      double phloem_capacity = parameters.percentage_phloem_storage * out.culm_growth.phloem[day-1];

      // --- STEP 2: Compute per-organ excess (above storage limits) ---
      double excess_needles = std::max(sugar.needles - needle_capacity, 0.0);
      double excess_phloem  = std::max(sugar.phloem  - phloem_capacity, 0.0);

      // --- STEP 3: Total excess available for export ---
      double total_excess_available = excess_needles + excess_phloem;

      // --- STEP 4: Apply export rate ---
      double export_rate = 0.5; // You can adjust this efficiency
      double export_target = export_rate * total_excess_available;

      // --- STEP 5: Export from needles first ---
      double export_from_needles = std::min(sugar.needles, std::min(excess_needles, export_target));
      sugar.needles -= export_from_needles;
      export_target -= export_from_needles;

      // --- STEP 6: Export remaining from phloem (if needed) ---
      double export_from_phloem = std::min(sugar.phloem, std::min(excess_phloem, export_target));
      sugar.phloem -= export_from_phloem;
      export_target -= export_from_phloem;

      // --- STEP 7: Finalize export ---
      double total_exported = export_from_needles + export_from_phloem;
      sugar_out_of_system += total_exported;
      sugar.surplus = total_exported;  // For tracking in output

      /*
       * Balence Check
       */

      double carbo_ending_3 = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
        starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + sugar_out_of_system;
      double difference_3 = carbo_beginning - carbo_ending_3;
      if (difference_3 > 1e-8) {
        std::cout << "On day " << day + 1 << " After Excess: More carbohydrates compared to the BEGINNING of the day. Difference is: " << difference_3 << "\n";
      } else if (difference_3 < -1e-8) {
        std::cout << "On day " << day + 1 << " After Excess: More carbohydrate compared to the END of the day. Difference is: " << difference_3 << "\n";
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

        // TODO: make these input parameters
        double N = 0.5;
        double mycorrhizal_nitrogen_demand = 0.2;

        // TODO: make this into an input parameter
        // TODO: consdier a different uptake rate for the sugar surplus and the sugar.mycorrhiza!
        uptake = nitrogen_uptake(N,
                                 sugar.mycorrhiza,
                                 sugar.surplus,
                                 out.culm_growth.mycorrhiza[day-1],
                                 out.culm_growth.roots[day-1],
                                 mycorrhizal_nitrogen_demand,
                                 TRUE);

        // C:N ratios are from the Korhonen 2013 paper
        // TODO: real value for the roots
        // The parameters here are the nitrogen content of each organ per biomass divided by 0.5 to get N:C

        // TODO: why does the height have more nitrogen than the wall? This doesn't seem quite right...

        nitrogen_balance = nitrogen_balance + uptake.total_uptake -
          4.9/0.5 * std::min(storage_term.needles, nitrogen_capacity.needles) * tree_state.needles -
          2.1/0.5 * std::min(storage_term.needles, nitrogen_capacity.bud) * tree_state.bud -
          3.7/0.5 * phloem_growth_share * (std::min(storage_term.phloem, nitrogen_capacity.wall) * tree_state.wall -
          4.9/0.5 * std::min(storage_term.phloem, nitrogen_capacity.height) * tree_state.height) -
          3.7/0.5 * xylem_sh_growth_share * (std::min(storage_term.xylem_st, nitrogen_capacity.wall) * tree_state.wall -
          4.9/0.5 * std::min(storage_term.xylem_st, nitrogen_capacity.height) * tree_state.height) -
          3.0/0.5 * std::min(storage_term.roots, nitrogen_capacity.roots) * tree_state.roots;
      }

      /*
       * Storage check
       */

      if ((starch.needles == 0) && (starch.phloem == 0) && (starch.xylem_sh == 0) && (starch.xylem_st == 0) && (starch.roots == 0) &&
          (sugar.needles == 0) && (sugar.phloem == 0) && (sugar.xylem_sh == 0) && (sugar.xylem_st == 0) && (sugar.roots == 0)) {
        std::cerr << " Day " << day << " No total storage - plant died!" << "\n";
        // tree_alive = false;
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
      // TODO: put this in initalisation as no longer going through on day 0
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

      if ((tree_state.sH > parameters.sHc) & ((sugar.needles + starch.needles) > 0.07)) {
        myco_allocation = PF * 0.3;
      } else {
        myco_allocation = 0.0;
      }

      sugar_all = sugar.needles + PF - tree_state.use + tree_state.release - std::min(storage_term.respiration, nitrogen_capacity_all) * tree_state.Rm_a -
        (1 + common.Rg_S) * std::min(storage, nitrogen_capacity_all) * tree_state.height -
        (1 + common.Rg_S) * std::min(storage, nitrogen_capacity_all) * tree_state.diameter -
        (1 + common.Rg_N) * std::min(storage, nitrogen_capacity_all) * tree_state.needles -
        (1 + common.Rg_R) * std::min(storage, nitrogen_capacity_all) * tree_state.roots -
        (1 + common.Rg_N) * std::min(storage, nitrogen_capacity_all) * tree_state.bud -
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
      tree_alive = false;
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

    respiration_growth = (common.Rg_S) * std::min(storage, nitrogen_capacity_all) * tree_state.height -
      (common.Rg_S) * std::min(storage, nitrogen_capacity_all) * tree_state.diameter -
      (common.Rg_N) * std::min(storage, nitrogen_capacity_all) * tree_state.needles -
      (common.Rg_R) * std::min(storage, nitrogen_capacity_all) * tree_state.roots -
      (common.Rg_N) * std::min(storage, nitrogen_capacity_all) * tree_state.bud;
    respiration_maintainence = storage_term.respiration * tree_state.Rm_a;
  }

  /*
   * OUTPUT!
   */

  log_sugar(day,
            days_gone,
            sugar,
            starch,
            storage_term,
            nitrogen_capacity,
            respiration_growth,
            respiration_maintainence,
            nitrogen_balance,
            uptake,
            tree_alive,
            out);
}
