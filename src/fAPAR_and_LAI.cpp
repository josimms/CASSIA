#include "CASSIA.h"

// The function definition
void compute_fAPAR_used(int day,
                        int days_gone,
                        double LAI,
                        double max_needles,
                        const Settings& boolsettings,
                        const CASSIA_parameters& parameters,
                        const growth_state& tree_state,
                        const photosynthesis_out& photosynthesis,
                        const weather_all& climate,
                        output_vector& all_out
) {
  double fAPAR_used = 0.0;
  double LAI_within_year = 0.0;

  if (!boolsettings.photosynthesis_as_input && boolsettings.fAPAR_Tian) {
    // Tian 2021 method
    double f_modifer = (max_needles == 0.0) ? 0.0 : tree_state.needles / max_needles;

    if (day < 182) {
      LAI_within_year = LAI * (parameters.n_age - 1) / parameters.n_age +
        (1.0 / parameters.n_age) * f_modifer * LAI;
    } else if (day > 244) {
      LAI_within_year = LAI * (parameters.n_age - 1) / parameters.n_age +
        (1.0 / parameters.n_age) * photosynthesis.fS * LAI;
    } else {
      LAI_within_year = LAI;
    }

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
    fAPAR_used = climate.fAPAR[day];
  } else {
    fAPAR_used = 0.0;
  }

  // Logging
  all_out.LAI[days_gone + day] = LAI_within_year;
  all_out.fAPAR[days_gone + day] = fAPAR_used;
}
