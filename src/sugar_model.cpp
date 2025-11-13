#include "CASSIA.h"

/*
 * == Check for NaN ==
 */

void check_tracker_for_nan(const carbo_tracker& tracker, const std::string& name, int day) {
  auto check_value = [&](double val, const std::string& label) {
    if (std::isnan(val) || std::isinf(val)) {
      std::cout << "⚠️ Day " << day << " - " << name << "." << label
                << " is NaN or Inf! Value: " << val << "\n";
    } else if (val == 0.0) {
      std::cout << "⚠️ Day " << day << " - " << name << "." << label
                << " is ZERO. Value: " << val << "\n";
    } else if (val < 0.0) {
      std::cout << "⚠️ Day " << day << " - " << name << "." << label
                << " is NEGATIVE! Value: " << val << "\n";
    }
  };

  check_value(tracker.needles,  " needles");
  check_value(tracker.phloem,   " phloem");
  check_value(tracker.roots,    " roots");
  check_value(tracker.xylem_sh, " xylem_sh");
  check_value(tracker.xylem_st, " xylem_st");
}


/*
 * === Sugar Storage ===
 */

double storage_update_organs(double growth, double sugar, double lower_limit, bool tree_alive) {
  if (!tree_alive || std::isnan(sugar)) return 0.0;

  if (growth <= 0.0) return 0.0; // avoid division by zero or negative growth

  // Sugar available beyond the minimum survival level
  double available = sugar - lower_limit;

  // fraction of growth that can be satisfied
  double fraction = 0.0;
  if (available > 0.0) {
    fraction = available / growth;
    fraction = std::min(fraction, 1.0);   // cap at 1
  }

  fraction = std::max(fraction, 0.0);       // ensure non-negative if small negative floating point

  return fraction;
}


/*
 * === Respiration sugar check, capped at a lower bound ===
 */

double respiration_check(double storage_term, double lower_bound) {
  double out = storage_term;
  if (storage_term < lower_bound) {
    out = 0;
  }
  return(out);
}

/*
 * === Nitrogen Storage ===
 */

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
    storage_capcity = 4.5;
  } else if (organ == "height") {
    storage_capcity = 4.5;
  } else if (organ == "roots") {
    storage_capcity = 3.0;            // TODO: find a better value for this!
  } else if (organ == "all") {
    // TODO: Average of the biomass
    storage_capcity = 4.0;
  }

  if (nitrogen_balance > storage_capcity) {
    out = 1.0;
  } else {
    out = std::min(nitrogen_balance/storage_capcity, 1.0);
  }

  return(out);
}

/*
 *  === To Starch ===
 */

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
  double zeta = 5;

  // Sugar surplus after winter needs
  double S = sugar - winter;
  double sugar_limitation = sigma(S, zeta);
  double starch_limitation = sigma(starch, zeta);
  double exponent = max_capacity - S - starch;
  double max_capacity_limitation = sigma(exponent, zeta);
  double synthesis = Ad0 * S * max_capacity_limitation;
  double remobilisation = delta * starch * starch_limitation;
  double proposed_flux = sugar_limitation * synthesis - (1 - sugar_limitation) * remobilisation;

  // Print all intermediate values in one line
  // std::cout
  // << "DEBUG safe_to_starch_transfer: "
  // << "sugar=" << sugar << " "
  // << "winter=" << winter << " "
  // << "S=" << S << " "
  // << "sugar_limitation=" << sugar_limitation << " "
  // << "starch=" << starch << " "
  // << "starch_limitation=" << starch_limitation << " "
  // << "max_capacity=" << max_capacity << " "
  // << "exponent=" << exponent << " "
  // << "max_capacity_limitation=" << max_capacity_limitation << " "
  // << "Ad0=" << Ad0 << " "
  // << "synthesis=" << synthesis << " "
  // << "delta=" << delta << " "
  // << "remobilisation=" << remobilisation << " "
  // << "proposed_flux=" << proposed_flux
  // << std::endl;

  return proposed_flux;
}


/*
 * === Filter for Nitrogen uptake to Roots ===
 */

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

/*
 * === Nitrogen Uptake ===
 */

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

/*
 * === Sugar transfer to mycorrhiza if the excess sugar transfer is worked out earlier ===
 */

double constant_excess(double N_T, double N_r, double N_m,
                       double excess, double beta_m) {
  // Quadratic coefficients
  double A = N_T - N_m;
  double B = N_T * (excess + beta_m + 1e-6) - excess * N_r - excess * N_m - N_m * 1e-6;
  double C = N_T * beta_m * (excess + 1e-6) - excess * N_r * beta_m;

  // Calculate discriminant
  double discriminant = B*B - 4*A*C;
  if (discriminant < 0) {
    std::cerr << "No real solution exists.\n";
    return 0.0;
  }

  double sqrt_disc = std::sqrt(discriminant);
  double x1 = (-B + sqrt_disc) / (2*A);
  double x2 = (-B - sqrt_disc) / (2*A);

  // Return only the positive solution
  if (x1 > 0 && x2 > 0) {
    std::cout << "Two positive solutions returned the minimum.";
    return std::min(x1, x2);
  }
  if (x1 > 0) return x1;
  if (x2 > 0) return x2;

  std::cerr << "No positive solution exists.\n";
  return 0.0;
}

/*
 * === Emergency sugar transfer ===
 */

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
 * === Calculate Movement ===
 */

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
 * ============================= SUGAR MODEL ===============================
 */

void sugar_model(int day,
                 int days_gone,

                 double TAir,
                 double PAR,
                 double PF,

                 CASSIA_parameters parameters,
                 CASSIA_common common,

                 bool nitrogen_change,
                 bool nitrogen_contrast,
                 bool mycorrhiza_passive,
                 bool surplus_c,
                 bool tree_alive,
                 bool new_sugar_model,



                 const growth_state& tree_state,
                 double& nitrogen_balance,
                 uptake_structre& uptake,
                 carbo_tracker& sugar,
                 carbo_tracker& starch,
                 carbo_tracker& storage_term,
                 growth_out& nitrogen_capacity,
                 output_vector& out) {

  double respiration_growth = 0.0;
  double respiration_maintainence = 0.0;

  if (new_sugar_model) {
    // check_tracker_for_nan(sugar, "Sugar, before anything", day + days_gone + 1);
    // check_tracker_for_nan(starch, "Starch, before anything", day + days_gone + 1);

    double denom = out.culm_growth.phloem[day + days_gone] + out.culm_growth.xylem_sh[day + days_gone] + out.culm_growth.xylem_st[day + days_gone];

    double phloem_share   = denom > 0 ? out.culm_growth.phloem[day + days_gone] / denom : 0.0;
    double xylem_st_share = denom > 0 ? out.culm_growth.xylem_st[day + days_gone] / denom : 0.0;
    double xylem_sh_share = denom > 0 ? out.culm_growth.xylem_sh[day + days_gone] / denom : 0.0;

    /*
     * Wintering
     */
    // 12% to 22% of water content reduces the freezing rate of water
    // Just using the values from autumn at the moment.
    // Temperature and PAR average used to control the decrease / increase in allocation.

    double winter = (1.0/(1.0 + std::exp(0.5*(TAir-5))) + 1.0/(1.0 + std::exp(0.2*(PAR-30))))/2.0;

    carbo_tracker winter_costs{};
    winter_costs.needles  = 0.0536 * winter + parameters.lower_bound_needles;
    winter_costs.phloem   = 0.2736 * winter + parameters.lower_bound_phloem;
    winter_costs.roots    = 0.0300 * winter + parameters.lower_bound_roots;
    winter_costs.xylem_st = 0.0569 * winter + parameters.lower_bound_xylem_st;
    winter_costs.xylem_sh = 0.0231 * winter + parameters.lower_bound_xylem_sh;

    /*
     * Storage
     */

    storage_term.needles  = storage_update_organs(tree_state.RmN + (1 + common.Rg_N) * (tree_state.needles + tree_state.bud), sugar.needles, parameters.lower_bound_needles, tree_alive);
    storage_term.phloem   = storage_update_organs(phloem_share * tree_state.RmS + phloem_share * (1 + common.Rg_S) * (tree_state.wall + tree_state.height) + tree_state.use - tree_state.release, sugar.phloem, parameters.lower_bound_phloem, tree_alive);
    storage_term.xylem_sh = storage_update_organs(xylem_sh_share * tree_state.RmS + xylem_sh_share * (1 + common.Rg_S) * (tree_state.wall + tree_state.height), sugar.xylem_sh, parameters.lower_bound_xylem_sh, tree_alive);
    storage_term.xylem_st = storage_update_organs(xylem_st_share * tree_state.RmS + xylem_st_share * (1 + common.Rg_S) * (tree_state.wall + tree_state.height), sugar.xylem_st, parameters.lower_bound_xylem_st, tree_alive);
    storage_term.roots    = storage_update_organs(tree_state.RmR + (1 + common.Rg_R) * tree_state.roots, sugar.roots, parameters.lower_bound_roots, tree_alive);

    // MAKING A BALANCE RELATIONSHIP

    // Note the relationship now is using the all component as the all reference for the N:C balance is a weight based average of the different components.
    // The storage is in different compartments, however so there could be an expansion here!
    nitrogen_capacity.needles = nitrogen_storage(nitrogen_balance, "all");
    nitrogen_capacity.bud     = nitrogen_storage(nitrogen_balance, "all");
    nitrogen_capacity.wall    = nitrogen_storage(nitrogen_balance, "all");
    nitrogen_capacity.height  = nitrogen_storage(nitrogen_balance, "all");
    nitrogen_capacity.roots   = nitrogen_storage(nitrogen_balance, "all");

    /*
     * Carbon allocation plan growth vs allocation
     *
     * Note the allocation is here as this means that the growth and allocation are worked out before the growth is applied.
     * Although the nitrogen uptake affect will be seen only in the next iteration thanks to the uptake of storage being calculated before this stage and nitrogen as the last stage in the algorithm
     */

    // A FIXED ALLOCATION ASSUMED
    double total_sugar = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots;
    double myco_allocation = parameters.mycorrhiza_threshold * total_sugar;

    if (surplus_c) {
      // STRATEGY 1: This strategy assumes that only the surplus is given to the mycorrhiza.
      myco_allocation = 0.0;
    } else if (nitrogen_contrast) {
      // Calculate the difference
      double growth_defficit = storage_term.needles - nitrogen_capacity.needles +
        storage_term.needles - nitrogen_capacity.bud +
        storage_term.roots - nitrogen_capacity.roots +
        phloem_share * (storage_term.phloem - nitrogen_capacity.wall) + xylem_sh_share * (storage_term.xylem_sh - nitrogen_capacity.wall) + xylem_st_share * (storage_term.xylem_st - nitrogen_capacity.wall) +
        phloem_share * (storage_term.phloem - nitrogen_capacity.height) + xylem_sh_share * (storage_term.xylem_sh - nitrogen_capacity.height) + xylem_st_share * (storage_term.xylem_st - nitrogen_capacity.height);

      if (growth_defficit <= 0.0) {
        // Less sugar, so sugar is limiting the growth
        myco_allocation = 0.0;
      } else if (growth_defficit > 0.0) {
        // More sugar so nitrogen is limiting growth

        // The missing nitrogen as the capacity difference is multiplied by the nitrogen needed for the growth
        double nitrogen_target = 4.9/0.5 * (storage_term.needles - nitrogen_capacity.needles) * tree_state.needles -
          2.1/0.5 * (storage_term.needles - nitrogen_capacity.bud) * tree_state.bud -
          4.5/0.5 * phloem_share * ((storage_term.phloem - nitrogen_capacity.wall) * tree_state.wall +
          (storage_term.phloem - nitrogen_capacity.height) * tree_state.height) -
          4.5/0.5 * xylem_sh_share * ((storage_term.xylem_st - nitrogen_capacity.wall) * tree_state.wall +
          (storage_term.xylem_sh - nitrogen_capacity.height) * tree_state.height) -
          4.5/0.5 * xylem_st_share * ((storage_term.xylem_st - nitrogen_capacity.wall) * tree_state.wall +
          (storage_term.xylem_st - nitrogen_capacity.height) * tree_state.height) -
          3.0/0.5 * (storage_term.roots - nitrogen_capacity.roots) * tree_state.roots;

        // The sugar that could have gone into growth is seen as the investment pool
        double sugar_to_invest = (1 + common.Rg_N) * ((storage_term.needles - nitrogen_capacity.needles) * tree_state.needles + (storage_term.needles - nitrogen_capacity.bud) * tree_state.bud) +
          (1 + common.Rg_R) * (storage_term.roots - nitrogen_capacity.roots) * tree_state.roots +
          phloem_share * (1 + common.Rg_S) * ((storage_term.phloem - nitrogen_capacity.wall) * tree_state.wall + (storage_term.phloem - nitrogen_capacity.height) * tree_state.height) +
          xylem_st_share * (1 + common.Rg_S) * ((storage_term.xylem_st - nitrogen_capacity.wall) * tree_state.wall + (storage_term.xylem_st - nitrogen_capacity.height) * tree_state.height) +
          xylem_sh_share * (1 + common.Rg_S) * ((storage_term.xylem_sh - nitrogen_capacity.wall) * tree_state.wall + (storage_term.xylem_sh - nitrogen_capacity.height) * tree_state.height);

        // Excess is consisted constant here and the same as the previous iteration.
        double half_sat_myco = 0.5;
        // TODO: this should be possible, because the uptake is given as an output to the model
        myco_allocation = std::min(constant_excess(nitrogen_target, uptake.root_uptake, uptake.ectomycorrhizal_uptake, sugar.surplus, half_sat_myco), sugar_to_invest);
      }
    }

    /*
     * Growth and respiration
     */

    // CALCULATE THE GROWTH WITH THE LAST ITERATIONS SUGAR AND NITROGEN
    carbo_tracker growth_resp{};
    growth_resp.needles =                   tree_state.RmN * storage_term.needles  +                  (1 + common.Rg_N) * (std::min(storage_term.needles, nitrogen_capacity.needles)  * tree_state.needles + std::min(storage_term.needles, nitrogen_capacity.bud) * tree_state.bud);
    growth_resp.roots =                     tree_state.RmR * storage_term.roots    +                  (1 + common.Rg_R) * std::min(storage_term.roots, nitrogen_capacity.roots)       * tree_state.roots;
    growth_resp.phloem =     phloem_share * tree_state.RmS * storage_term.phloem   + phloem_share   * (1 + common.Rg_S) * (std::min(storage_term.phloem, nitrogen_capacity.wall)      * tree_state.wall    + std::min(storage_term.phloem, nitrogen_capacity.height) * tree_state.height)   + std::min(storage_term.phloem, nitrogen_capacity.height) * (tree_state.use - tree_state.release);
    growth_resp.xylem_st = xylem_st_share * tree_state.RmS * storage_term.xylem_st + xylem_st_share * (1 + common.Rg_S) * (std::min(storage_term.xylem_st, nitrogen_capacity.wall)    * tree_state.wall    + std::min(storage_term.xylem_st, nitrogen_capacity.height) * tree_state.height);
    growth_resp.xylem_sh = xylem_sh_share * tree_state.RmS * storage_term.xylem_sh + xylem_sh_share * (1 + common.Rg_S) * (std::min(storage_term.xylem_sh, nitrogen_capacity.wall)    * tree_state.wall    + std::min(storage_term.xylem_sh, nitrogen_capacity.height) * tree_state.height);

    // APPLY THE GROWTH
    double before_growth = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots + PF;

    // TODO: problem, this is the last iteration of the sugar and the new iteration in the nitrogen
    if (sugar.needles < growth_resp.needles) {
      std::cout << "Day " << day + days_gone + 1 << ": Stoage term is not capping growth enough needles. Needles storage " << storage_term.needles << " sugar.needles " << sugar.needles << " growth_resp.needles " << growth_resp.needles << "\n";
    }
    sugar.needles = sugar.needles + PF - growth_resp.needles;

    if (sugar.phloem < growth_resp.phloem) {
      std::cout << "Day " << day + days_gone + 1 << ": Stoage term is not capping growth enough phloem. Phloem storage " << storage_term.phloem << " sugar.phloem " << sugar.phloem << " growth_resp.phloem " << growth_resp.phloem << "\n";
    }
    sugar.phloem = sugar.phloem - growth_resp.phloem;


    if (sugar.roots < growth_resp.roots) {
      std::cout << "Day " << day + days_gone + 1 << ": Stoage term is not capping growth enough roots. Root storage " << storage_term.roots << " sugar.roots " << sugar.roots << " growth_resp.roots " << growth_resp.roots << "\n";
    }
    sugar.roots = sugar.roots - growth_resp.roots;


    if (sugar.xylem_sh < growth_resp.xylem_sh) {
      std::cout << "Day " << day + days_gone + 1 << ": Stoage term is not capping growth enough xylem sh. Xylem Sh storage " << storage_term.xylem_sh  << " sugar.xylem_sh " << sugar.xylem_sh << " growth_resp.xylem_sh " << growth_resp.xylem_sh << "\n";
    }
    sugar.xylem_sh = sugar.xylem_sh - growth_resp.xylem_sh;


    if (sugar.xylem_st < growth_resp.xylem_st) {
      std::cout << "Day " << day + days_gone + 1 << ": Stoage term is not capping growth enough xylem st. Xylem St storage " << storage_term.xylem_st  << " sugar.xylem_st " << sugar.xylem_st << " growth_resp.xylem_st " << growth_resp.xylem_st << "\n";
    }
    sugar.xylem_st = sugar.xylem_st - growth_resp.xylem_st;

    // CHECKS
    double after_growth = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
      growth_resp.needles + growth_resp.roots + growth_resp.phloem + growth_resp.xylem_st + growth_resp.xylem_sh;

    if (after_growth - before_growth > 1e-7) {
      std::cout << "Warning: Growth pools are not consistant, after growth is more than before!\n";
    } else if (after_growth - before_growth < -1e-7) {
      std::cout << "Warning: Growth pools are not consistant, beofre growth is more than after!\n";
    }

    // check_tracker_for_nan(sugar, "Sugar, after growth", day + days_gone + 1);
    // check_tracker_for_nan(starch, "Starch, after growth", day + days_gone + 1);

    // NITROGEN BALANCE CHANGE IF LINKED WITH GROWTH

    if (nitrogen_change) {
      // UPDATE THE NITROGEN BALANCE
      // C:N ratios are from the Korhonen 2013 paper
      // TODO: real value for the roots
      // The parameters here are the nitrogen content of each organ per biomass divided by 0.5 to get N:C
      nitrogen_balance = nitrogen_balance + uptake.total_uptake -
        4.9/0.5 * std::min(storage_term.needles, nitrogen_capacity.needles) * tree_state.needles -
        2.1/0.5 * std::min(storage_term.needles, nitrogen_capacity.bud) * tree_state.bud -
        4.5/0.5 * phloem_share * (std::min(storage_term.phloem, nitrogen_capacity.wall) * tree_state.wall +
        std::min(storage_term.phloem, nitrogen_capacity.height) * tree_state.height) -
        4.5/0.5 * xylem_st_share * (std::min(storage_term.xylem_st, nitrogen_capacity.wall) * tree_state.wall +
        std::min(storage_term.xylem_st, nitrogen_capacity.height) * tree_state.height) -
        4.5/0.5 * xylem_sh_share * (std::min(storage_term.xylem_sh, nitrogen_capacity.wall) * tree_state.wall +
        std::min(storage_term.xylem_sh, nitrogen_capacity.height) * tree_state.height) -
        3.0/0.5 * std::min(storage_term.roots, nitrogen_capacity.roots) * tree_state.roots;

      nitrogen_capacity.height = nitrogen_balance;
    }

    /*
     * Allocation based on growth sink demand
     */
    // Check
    double before_transfer = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
       starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots;

    // --- INPUTS ---
    double sugar_needles = sugar.needles;
    double sugar_phloem  = sugar.phloem;

    // --- STEP 1: Compute sink demands and capacities ---

    // ROOTS
    double root_storage_total = sugar.roots + starch.roots;
    double root_capacity_max = parameters.percentage_roots_storage * out.culm_growth.roots[day + days_gone - 1];
    double root_capacity_remaining = std::max(root_capacity_max - root_storage_total, 0.0);
    double demand_roots = growth_resp.roots + winter_costs.roots;
    double actual_root_demand = std::min(demand_roots, root_capacity_remaining);

    // XYLEM ST
    double xylem_st_storage_total = sugar.xylem_st + starch.xylem_st;
    double xylem_st_capacity_max = parameters.percentage_xylem_st_storage * out.culm_growth.xylem_st[day + days_gone];
    double xylem_st_capacity_remaining = std::max(xylem_st_capacity_max - xylem_st_storage_total, 0.0);
    double demand_xylem_st = growth_resp.xylem_st + winter_costs.xylem_st;
    double actual_xylem_st_demand = std::min(demand_xylem_st, xylem_st_capacity_remaining);

    // XYLEM SH
    double xylem_sh_storage_total = sugar.xylem_sh + starch.xylem_sh;
    double xylem_sh_capacity_max = parameters.percentage_xylem_sh_storage * out.culm_growth.xylem_sh[day + days_gone];
    double xylem_sh_capacity_remaining = std::max(xylem_sh_capacity_max - xylem_sh_storage_total, 0.0);
    double demand_xylem_sh = growth_resp.xylem_sh + winter_costs.xylem_sh;
    double actual_xylem_sh_demand = std::min(demand_xylem_sh, xylem_sh_capacity_remaining);

    // PHLOEM SELF-DEMAND
    double demand_phloem = growth_resp.phloem + winter_costs.phloem;

    // Compute total effective (capacity-limited + flow-through + phloem) demand
    double capacity_limited_demand = actual_root_demand + actual_xylem_st_demand + actual_xylem_sh_demand;
    double total_effective_demand = capacity_limited_demand + myco_allocation + demand_phloem;

    // --- STEP 2: Compute phloem gap to meet demands ---
    // Note: winter costs are not counted here as they are part of the demand
    double available_phloem = sugar_phloem;
    double phloem_to_cover_sinks = std::max(total_effective_demand - available_phloem, 0.0);

    // --- STEP 3: Determine needle-to-phloem transfer (capped by reserve and phloem need) ---
    // Note: 1.2 added to give a buffer
    double needle_reserve = std::max(sugar_needles - 2.0 * winter_costs.needles - growth_resp.needles, 0.0);
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
    double myco_export = std::min(phloem_available_for_transfer, myco_allocation);
    sugar.phloem -= myco_export;
    phloem_available_for_transfer -= myco_export;

    // concentration_gradient.roots_to_myco = myco_export;
    sugar.mycorrhiza = myco_export;

    // --- BALANCE CHECK ---
    double after_transfer = sugar.needles + sugar.phloem + sugar.xylem_sh + sugar.xylem_st + sugar.roots +
    starch.needles + starch.phloem + starch.xylem_sh + starch.xylem_st + starch.roots + myco_export;

    double difference_2 = before_transfer - after_transfer;

    if (difference_2 > 1e-8) {
     std::cout << "Day " << day + days_gone + 1 << ": More carbs before internal transfer. Difference: " << difference_2 << "\n";
    } else if (difference_2 < -1e-8) {
     std::cout << "Day " << day + days_gone + 1 << ": More carbs after internal transfer: " << difference_2 << "\n";
    }

    // check_tracker_for_nan(sugar, "Sugar, after transfer", day + days_gone + 1);
    // check_tracker_for_nan(starch, "Starch, after transfer", day + days_gone + 1);

    /*
     * Starch
     */
    carbo_tracker to_starch{};

    // === NEEDLES ===
    {
      // TODO: leaf mass rather than needles
      double capacity = parameters.percentage_needle_storage * out.culm_growth.leaf_mass[day + days_gone];

      to_starch.needles = safe_to_starch_transfer(
        sugar.needles, winter_costs.needles, starch.needles,
        parameters.Ad0_needles, parameters.delta_needles, capacity
      );

      sugar.needles -= to_starch.needles;
      starch.needles += to_starch.needles;

      if (sugar.needles < 0.0) std::cout << "Day " << day + days_gone + 1 << " sugar.needles is negative after starch. Sugar " << sugar.needles << " to_starch.needles " << to_starch.needles << "\n";
      if (starch.needles < 0.0) std::cout << "Day " << day + days_gone + 1 << " starch.needles is negative after starch. Starch " << starch.needles << " to_starch.needles " << to_starch.needles << "\n";
    }

    // === PHLOEM ===
    {
      double capacity = parameters.percentage_phloem_storage * out.culm_growth.phloem[day + days_gone];

      to_starch.phloem = safe_to_starch_transfer(
        sugar.phloem, winter_costs.phloem, starch.phloem,
        parameters.Ad0_phloem, parameters.delta_phloem, capacity
      );

      sugar.phloem -= to_starch.phloem;
      starch.phloem += to_starch.phloem;

      if (sugar.phloem < 0.0) std::cout << "Day " << day + days_gone + 1 << " sugar.phloem is negative after starch. Sugar " << sugar.phloem << " to_starch.phloem " << to_starch.phloem << "\n";
      if (starch.phloem < 0.0) std::cout << "Day " << day + days_gone + 1 << " starch.phloem is negative after starch. Starch " << starch.phloem << " to_starch.phloem " << to_starch.phloem << "\n";
    }

    // === ROOTS ===
    {
      double capacity = parameters.percentage_roots_storage * out.culm_growth.roots[day + days_gone];

      to_starch.roots = safe_to_starch_transfer(
       sugar.roots, winter_costs.roots, starch.roots,
       parameters.Ad0_roots, parameters.delta_roots, capacity
      );

      sugar.roots -= to_starch.roots;
      starch.roots += to_starch.roots;

      if (sugar.roots < 0.0) std::cout << "Day " << day + days_gone + 1 << " sugar.roots is negative after starch. Sugar " << sugar.roots << " to_starch.roots " << to_starch.roots << "\n";
      if (starch.roots < 0.0) std::cout << "Day " << day + days_gone + 1 << " starch.roots is negative after starch. Starch " << starch.roots << " to_starch.roots " << to_starch.roots << "\n";
    }

    // === XYLEM SH ===
    {
      double capacity = parameters.percentage_xylem_sh_storage * out.culm_growth.xylem_sh[day + days_gone];

      to_starch.xylem_sh = safe_to_starch_transfer(
       sugar.xylem_sh, winter_costs.xylem_sh, starch.xylem_sh,
       parameters.Ad0_xylem_sh, parameters.delta_xylem_sh, capacity
      );

      sugar.xylem_sh -= to_starch.xylem_sh;
      starch.xylem_sh += to_starch.xylem_sh;

      if (sugar.xylem_sh < 0.0) std::cout << "Day " << day + days_gone + 1 << " sugar.xylem_sh is negative after starch. Sugar " << sugar.xylem_sh << " to_starch.xylem_sh " << to_starch.xylem_sh << "\n";
      if (starch.xylem_sh < 0.0) std::cout << "Day " << day + days_gone + 1 << " starch.xylem_sh is negative after starch. Starch " << starch.xylem_sh << " to_starch.xylem_sh " << to_starch.xylem_sh << "\n";
    }

    // === XYLEM ST ===
    {
      double capacity = parameters.percentage_xylem_st_storage * out.culm_growth.xylem_st[day + days_gone];

      to_starch.xylem_st = safe_to_starch_transfer(
       sugar.xylem_st, winter_costs.xylem_st, starch.xylem_st,
       parameters.Ad0_xylem_st, parameters.delta_xylem_st, capacity
      );

      sugar.xylem_st -= to_starch.xylem_st;
      starch.xylem_st += to_starch.xylem_st;

      if (sugar.xylem_st < 0.0) std::cout << "Day " << day + days_gone + 1 << " sugar.xylem_st is negative after starch. Sugar " << sugar.xylem_st << " to_starch.xylem_st " << to_starch.xylem_st << "\n";
      if (starch.xylem_st < 0.0) std::cout << "Day " << day + days_gone + 1 << " starch.xylem_st is negative after starch. Starch " << starch.xylem_st << " to_starch.xylem_st " << to_starch.xylem_st << "\n";
    }

    // check_tracker_for_nan(sugar, "Sugar, after starch", day + days_gone + 1);
    // check_tracker_for_nan(starch, "Starch, after starch", day + days_gone + 1);

    /*
     * === Excess sugar export beyond capacity ===
     */

    // --- STEP 1: Compute per-organ storage capacities ---
    double needle_capacity = parameters.percentage_needle_storage * out.culm_growth.needles[day + days_gone];
    double phloem_capacity = parameters.percentage_phloem_storage * out.culm_growth.phloem[day + days_gone];

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
    sugar.surplus = total_exported;  // For tracking in output

    // check_tracker_for_nan(sugar, "Sugar, after export", day + days_gone + 1);
    // check_tracker_for_nan(starch, "Starch, after export", day + days_gone + 1);

    /*
     * Uptake
     */
    if (nitrogen_change) {
      // UPTAKAE: TODO: root_biomass and mycorrhizal_biomass from last iteration

      double N = 0.5;
      double mycorrhizal_nitrogen_demand = 0.2;

      uptake = nitrogen_uptake(N,
                               myco_allocation,
                               sugar.surplus,
                               out.culm_growth.mycorrhiza[day + days_gone -1],
                               out.culm_growth.roots[day + days_gone -1],
                               mycorrhizal_nitrogen_demand,
                               mycorrhiza_passive);

    }
  } else {
    double respiration_growth = 0.0;
    double respiration_maintainence = 0.0;

    // Model
    double ak = 1 / (1 - 1/exp(parameters.alfa * (parameters.sugar00 + parameters.starch00 - parameters.Wala)));
    double storage{}, storage_term_Rm{}, sugar_all{}, starch_all{}, to_sugar{}, to_starch{};
    double nitrogen_capacity_all{};
    double myco_allocation{};
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
