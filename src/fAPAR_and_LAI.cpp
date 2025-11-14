#include "CASSIA.h"

// The function definition
void compute_fAPAR_used(int day,
                        int days_gone,
                        double LAI,
                        double max_needles,
                        bool& fS_reached_one,
                        bool& fN_reached_one,
                        repola_out repola,
                        const Settings& boolsettings,
                        const CASSIA_parameters& parameters,
                        const growth_state& tree_state,
                        const photosynthesis_out& photosynthesis,
                        const weather_all& climate,
                        output_vector& all_out
) {

  /*
   * Indexes
   */

  int index_ref = days_gone + day - 1;
  if (index_ref < 0) {
    index_ref = 0;
  }

  /*
   * Parameters
   */

  double fAPAR_used = 0.0;
  double LAI_within_year = 0.0;
  double LMA = 0.06; // 0.086; // kg / m-2
  double leaf_to_sapwood = 1.0 / 1400.0; // m2 m-2
  double needle_mass = repola.needle_mass; // TODO units?
  double needle_mass_sapwood = repola.needle_mass;

  // Logging
  all_out.LAI[days_gone + day] = LAI_within_year;
  all_out.fAPAR[days_gone + day] = fAPAR_used;

  // Tian 2021 method
  double culm_needle = 0.0;
  if (day > 0) {
    culm_needle = all_out.culm_growth.needles[index_ref];
  }
  double f_modifer = culm_needle/max_needles;
  if (f_modifer > 0.99) fN_reached_one = true;
  if (fN_reached_one) f_modifer = 1.0;

  double senescence = 0.0;
  // Day 50 added below as the first few days can give 1 which ruins the logic!
  if (photosynthesis.fS >= 1.0 && day > 50) fS_reached_one = true;
  if (fS_reached_one) senescence = photosynthesis.fS;

  // Calculated first as needle_mass changed otherwise
  needle_mass_sapwood = needle_mass * (parameters.n_age - 1.0) / parameters.n_age + (1.0 / parameters.n_age) * needle_mass * f_modifer;
  needle_mass         = needle_mass * (parameters.n_age - 1.0) / parameters.n_age + (1.0 / parameters.n_age) * needle_mass * f_modifer * senescence;

  if (!boolsettings.photosynthesis_as_input && boolsettings.fAPAR_Tian) {
    LAI_within_year = LAI * (parameters.n_age - 1.0) / parameters.n_age + (1.0 / parameters.n_age) * f_modifer * LAI - (1.0 / parameters.n_age) * senescence * LAI;

    fAPAR_used = 1.0 - std::exp(-0.52 * LAI_within_year);

    if (std::isnan(LAI_within_year)) {
      std::cout << "Day " << day << " LAI_within_year was NaN so replaced. TODO: fix\n";
      fAPAR_used = 0.7;
    }

    if (std::isnan(fAPAR_used)) {
      std::cout << "Day " << day << " fAPAR_used was NaN so replaced. TODO: fix\n";
      fAPAR_used = 0.7;
    }

  } else if (!boolsettings.photosynthesis_as_input && !boolsettings.fAPAR_Tian && boolsettings.preles) {
    fAPAR_used = climate.fAPAR[day + days_gone];
  } else {
    fAPAR_used = 0.0;
  }

  /*
   * Logging
   */

  // Leaf

  all_out.LAI[days_gone + day] = LAI_within_year;
  all_out.fAPAR[days_gone + day] = fAPAR_used;

  all_out.culm_growth.leaf_mass[day + days_gone] = needle_mass;
  all_out.culm_growth.leaf_area[day + days_gone] = needle_mass_sapwood / LMA;

  double sapwood_area = leaf_to_sapwood * all_out.culm_growth.leaf_area[day + days_gone];

  // Wood growth

  // (Scheistl Aalto, 2019): Mean wood density 200 kg C m−3
  all_out.culm_growth.sapwood[day + days_gone] = std::max(sapwood_area * all_out.culm_growth.height[index_ref] * 200.0, all_out.culm_growth.sapwood[index_ref]);

  // (Scheistl Aalto, 2019): "Sapwood was further divided to
  // 1) smaller branches and 2) bigger branches and trunk with ratio 1/9"
  all_out.culm_growth.xylem_sh[day + days_gone] = (1.0 / 9.0) * all_out.culm_growth.sapwood[day + days_gone];

  all_out.culm_growth.xylem_st[day + days_gone] = (8.0 / 9.0) * all_out.culm_growth.sapwood[day + days_gone];

  all_out.culm_growth.phloem[day + days_gone] = 0.1 * all_out.culm_growth.sapwood[day + days_gone];
}
