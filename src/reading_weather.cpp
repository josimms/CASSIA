#include "CASSIA.h"

weather_all readWeatherVariables(const Rcpp::DataFrame& weather, bool spp, bool preles, bool phydro) {
  weather_all weatherVariables;

  /*
   * Read the data
   */

  // Extract variables from the DataFrame
  weatherVariables.TAir = Rcpp::as<std::vector<double>>(weather["T"]);
  weatherVariables.TSoil_A = Rcpp::as<std::vector<double>>(weather["TSA"]);
  weatherVariables.TSoil_B = Rcpp::as<std::vector<double>>(weather["TSB"]);
  weatherVariables.Soil_Moisture = Rcpp::as<std::vector<double>>(weather["MB"]);
  weatherVariables.Precip = Rcpp::as<std::vector<double>>(weather["Rain"]);

  if (spp) {
    weatherVariables.Photosynthesis_IN = Rcpp::as<std::vector<double>>(weather["P"]);
  }

  if (preles) {
    weatherVariables.PAR = Rcpp::as<std::vector<double>>(weather["PAR"]);
    weatherVariables.VPD = Rcpp::as<std::vector<double>>(weather["VPD"]);
    weatherVariables.CO2 = Rcpp::as<std::vector<double>>(weather["CO2"]);
    weatherVariables.fAPAR = Rcpp::as<std::vector<double>>(weather["fAPAR"]);
  }

  if (phydro) {
    weatherVariables.PAR = Rcpp::as<std::vector<double>>(weather["PAR"]);
    weatherVariables.PAR_max = Rcpp::as<std::vector<double>>(weather["PAR_max"]);
    weatherVariables.VPD = Rcpp::as<std::vector<double>>(weather["VPD"]);
    weatherVariables.CO2 = Rcpp::as<std::vector<double>>(weather["CO2"]);
    weatherVariables.Nitrogen = Rcpp::as<std::vector<double>>(weather["Nitrogen"]);
    weatherVariables.PA = Rcpp::as<std::vector<double>>(weather["PA"]);
    weatherVariables.SWP = Rcpp::as<std::vector<double>>(weather["SWP"]);
  }

  /*
   * Print the weather data
   */

  // Print the first three values of each vector
  std::cout << "First three values of weather variables:" << std::endl;

  auto printFirstThree = [](const std::string& name, const std::vector<double>& values) {
    std::cout << name << ": ";
    for (size_t i = 0; i < std::min(values.size(), size_t(3)); ++i) {
      std::cout << values[i] << (i < 2 ? ", " : "\n");
    }
  };

  printFirstThree("TAir", weatherVariables.TAir);
  printFirstThree("TSoil_A", weatherVariables.TSoil_A);
  printFirstThree("TSoil_B", weatherVariables.TSoil_B);
  printFirstThree("Soil_Moisture", weatherVariables.Soil_Moisture);
  printFirstThree("Precip", weatherVariables.Precip);

  if (spp) {
    printFirstThree("Photosynthesis_IN", weatherVariables.Photosynthesis_IN);
  }

  if (preles || phydro) {
    printFirstThree("PAR", weatherVariables.PAR);
    printFirstThree("VPD", weatherVariables.VPD);
    printFirstThree("CO2", weatherVariables.CO2);
    printFirstThree("fAPAR", weatherVariables.fAPAR);
  }

  if (phydro) {
    printFirstThree("PAR_max", weatherVariables.PAR_max);
    printFirstThree("Nitrogen", weatherVariables.Nitrogen);
    printFirstThree("PA", weatherVariables.PA);
    printFirstThree("SWP", weatherVariables.SWP);
  }

  return weatherVariables;
}
