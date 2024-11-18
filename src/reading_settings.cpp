#include "CASSIA.h"
using namespace Rcpp;

// Function to convert an Rcpp::List to a Settings struct
Settings parseSettings(Rcpp::List settingsList) {
  Settings settings;

  // Access elements by name and assign to struct
  settings.storage_reset = Rcpp::as<bool>(settingsList["storage_reset"]);
  settings.storage_grows = Rcpp::as<bool>(settingsList["storage_grows"]);

  settings.LN_estim = Rcpp::as<bool>(settingsList["LN_estim"]);
  settings.mN_varies = Rcpp::as<bool>(settingsList["mN_varies"]);

  settings.LD_estim = Rcpp::as<bool>(settingsList["LD_estim"]);
  settings.sD_estim_T_count = Rcpp::as<bool>(settingsList["sD_estim_T_count"]);

  settings.LH_estim = Rcpp::as<bool>(settingsList["LH_estim"]);
  settings.trees_grow = Rcpp::as<bool>(settingsList["trees_grow"]);
  settings.growth_decreases = Rcpp::as<bool>(settingsList["growth_decreases"]);
  settings.needle_mass_grows = Rcpp::as<bool>(settingsList["needle_mass_grows"]);

  settings.phloem_trigger = Rcpp::as<bool>(settingsList["phloem_trigger"]);
  settings.mycorrhiza = Rcpp::as<bool>(settingsList["mycorrhiza"]);
  settings.root_as_Ding = Rcpp::as<bool>(settingsList["root_as_Ding"]);

  settings.sperling_model = Rcpp::as<bool>(settingsList["sperling_model"]);
  settings.myco_model = Rcpp::as<bool>(settingsList["myco_model"]);
  settings.xylogensis_option = Rcpp::as<bool>(settingsList["xylogenesis"]);

  settings.PRELES_GPP = Rcpp::as<bool>(settingsList["PRELES_GPP"]);
  settings.environmental_effect_xylogenesis = Rcpp::as<bool>(settingsList["environment_effect_xylogenesis"]);

  settings.photosynthesis_as_input = Rcpp::as<bool>(settingsList["photosynthesis_as_input"]);
  settings.phydro = Rcpp::as<bool>(settingsList["ecoevolutionary"]);

  settings.photoparameters = Rcpp::as<int>(settingsList["photoparameters"]);
  settings.temp_rise = Rcpp::as<bool>(settingsList["temp_rise"]);
  settings.drought = Rcpp::as<bool>(settingsList["drought"]);
  settings.Rm_acclimation = Rcpp::as<bool>(settingsList["Rm_acclimation"]);

  settings.CASSIA_graphs = Rcpp::as<bool>(settingsList["CASSIA_graphs"]);
  settings.tests = Rcpp::as<bool>(settingsList["tests"]);
  settings.etmodel = Rcpp::as<bool>(settingsList["etmodel"]);
  settings.LOGFLAG = Rcpp::as<bool>(settingsList["LOGFLAG"]);

  return settings;
}
