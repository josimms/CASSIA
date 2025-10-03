#include "CASSIA.h"

/*
 * === Sugar Storage ===
 */

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

    // Final result: limited to 1, capped so usage won't exceed available sugar
    out = std::min(1.0, ratio);
  } else {
    out = 0.0;
  }

  return out;
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
    return NAN;
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
  return NAN;
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

/*
 * Some of the parameters are site dependent, could define the site side in the main function code so I just import the right variables into this code
 *
 * needles mass per year, should define the year before this function!
 */

void sugar_model() {
  uptake_structre uptake;

  /*
   * Wintering
   */
  // 12% to 22% of water content reduces the freezing rate of water
  // Just using the values from autumn at the moment.
  // Temperature and PAR avaerage used to control the decrease / increase in allocation.

  double winter = (1.0/(1.0 + std::exp(0.5*(TAir-5))) + 1.0/(1.0 + std::exp(0.2*(PAR-30))))/2.0;

  carbo_tracker winter_costs;
  winter_costs.needles  = 0.0536 * winter + parameters.lower_bound_needles;
  winter_costs.phloem   = 0.2736 * winter + parameters.lower_bound_phloem;
  winter_costs.roots    = 0.0300 * winter + parameters.lower_bound_roots;
  winter_costs.xylem_st = 0.0569 * winter + parameters.lower_bound_xylem_st;
  winter_costs.xylem_sh = 0.0231 * winter + parameters.lower_bound_xylem_sh;

  /*
   * Uptake
   */

  /*
   * Storage
   */

  /*
   * Growth
   */

  /*
   * Allocation with storage
   */

  /*
   * Starch
   */

}
