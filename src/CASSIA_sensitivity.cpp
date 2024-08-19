#include "CASSIA.h"

std::vector<double> createRange(double start, double end, double step)
{
  std::vector<double> result;
  for (double i = start; i < end; i += step) {
    result.push_back(i);
  }
  return result;
};

// [[Rcpp::export]]
Rcpp::DataFrame replace_value_DataFrame(Rcpp::DataFrame df,
                                        double value,
                                        int ref)
{
  int nCols = df.cols();
  int nRows = df.rows();
  Rcpp::List result(nCols * nRows);
  result.attr("dim") = Rcpp::Dimension(nRows, nCols);
  colnames(result) = Rcpp::as<Rcpp::CharacterVector>(df.names());

  for (int i = 0; i < nCols; ++i) {
    std::vector<double> column = df[i];
    for (int j = 0; j < nRows; ++j) {
      double tmp = column[j];
      if (i == 0 && j == ref) {
          tmp = value;
        }
      result[i*nRows + j] = tmp;
    }
  }

  Rcpp::DataFrame df1(result);
  return df1;
}


/*
// [[Rcpp::export]]
Rcpp::List CASSIA_sensitivity(Rcpp::DataFrame bounds,
                              std::vector<std::string> names,

                                   int start_year,
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

                                   bool storage_rest,
                                   bool storage_grows,

                                   bool LH_estim,
                                   bool LN_estim,
                                   bool mN_varies,
                                   bool LD_estim,
                                   bool sD_estim_T_count,

                                   bool trees_grow,
                                   bool growth_decreases,
                                   bool needle_mass_grows,

                                   bool mycorrhiza,
                                   bool root_as_Ding,
                                   bool sperling_sugar_model,

                                   bool xylogensis_option,

                                   bool environmental_effect_xylogenesis,
                                   bool temp_rise,
                                   bool drought,
                                   bool Rm_acclimation,

                                   bool using_spp_photosynthesis,
                                   bool CASSIA_graphs,

                                   int etmodel,
                                   int LOGFLAG)
{
  // Create the vales, then plot in R
  Rcpp::List output;

  // Create bounds
  std::vector<double> lower_bounds = bounds["LL"];
  std::vector<double> upper_bounds = bounds["UL"];

  std::vector<std::string> parameters_names = pCASSIA_parameters.names();
  std::vector<std::string> sperling_names = pCASSIA_sperling.names();


  // Test bounds for a value
  for (double i = 0; i < lower_bounds.size(); i++) {
    // Create a range of values from the values in the bounds dataframe
    double step = (upper_bounds[i] - lower_bounds[i])/20;
    std::vector<double> range = createRange(lower_bounds[i], upper_bounds[i], step);

    Rcpp::List output_1;
    for (double value : range) {
      Rcpp::List out;
      std::vector<double> params_range = {22, 23, 24, 25, 26};
      if (std::find(params_range.begin(), params_range.end(), i) != params_range.end()) {
        // R references 62:66 for the parameters
        int ref = i + 39;
        Rcpp::DataFrame params = replace_value_DataFrame(pCASSIA_parameters, value, ref);
        // std::cout << parameters_names[ref] << " " << names[i] << "\n";
      } else {
        // R references c(50:54, 35:39, 40:44, 45:49,) for the parameters 1:22
        int ref;
        if (i < 5) {
          ref = i + 49;
        } else if (i >= 5 && i < 10) { // correct
          ref = i + 29;
        } else if (i >= 10 && i < 20) { // correct
          ref = i + 29;
        } else if (i == 20) {
          ref = i + 5;
        } else if (i == 21) {
          ref = i + 3;
        } else {
         std::cout << "ref outside range\n";
        }
        Rcpp::DataFrame sperling = replace_value_DataFrame(pCASSIA_sperling, value, ref);
        // std::cout << sperling_names[ref] << " " << names[i] << "\n";
      }

      out = CASSIA_yearly(start_year, end_year,
                          weather, GPP_ref,
                          pPREL, pCASSIA_parameters, pCASSIA_common, pCASSIA_ratios, pCASSIA_sperling, site,
                          needle_mass_in, Throughfall, storage_rest, storage_grows,
                          LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                          trees_grow, growth_decreases, needle_mass_grows,
                          mycorrhiza, root_as_Ding, sperling_sugar_model,
                          xylogensis_option, environmental_effect_xylogenesis,
                          temp_rise, drought, Rm_acclimation,
                          using_spp_photosynthesis, CASSIA_graphs,
                          etmodel, LOGFLAG);

      out.push_back(range);
      // Dynamically build list output
      output_1.push_back(out);
    }
    output.push_back(output_1);
  }

  return(output);
}
 */
