#include "CASSIA.h"
#include <cmath>

// **
// ** Gross and Net Assimilation
// **

// Hyytiala (the "Hyde" site) latitude, degrees North. The Plant-FATE/TreeLife
// companion codebase approximates this as 60.0 in its boreal config; CASSIA
// uses the real site value since it's calibrated specifically for Hyytiala.
static const double HYDE_LATITUDE_DEG = 61.85;

// Fraction of a 24h day that has daylight, from latitude and day-of-year.
// Ported 2026-10-29 from Assimilator::day_length_fraction() in the
// Plant-FATE/TreeLife companion codebase (src/assimilation.cpp), simplified
// to take an integer day-of-year directly (1 = Jan 1) instead of a decimal
// year, since CASSIA's day loop already has this as an integer. Replaces a
// hardcoded f_day_length = 0.5 that ignored boreal latitudes' huge seasonal
// swing (Hyytiala: ~6h in winter to ~19h in summer) - the growing season,
// when GPP matters most, was being systematically understated by that
// constant.
double day_length_fraction(double lat_deg, int day_of_year) {
  const double pi = M_PI;
  double lat = lat_deg * pi / 180.0;
  double z0 = 90.833 * pi / 180.0; // Zenith angle (sunrise/sunset)

  int N = day_of_year;

  // Solar declination (radians)
  double delta = -asin(0.39779 * cos(0.98565 * pi / 180.0 * (N + 10) +
                       1.914 * pi / 180.0 * sin(0.98565 * pi / 180.0 * (N - 2))));

  // Cosine of hour angle at sunrise/sunset
  double cos_h0 = (cos(z0) - sin(lat) * sin(delta)) / (cos(lat) * cos(delta));

  double D; // day length in hours
  if (cos_h0 <= -1.0) {
    D = 24.0; // Sun never sets
  } else if (cos_h0 >= 1.0) {
    D = 0.0; // Sun never rises
  } else {
    double H0 = acos(cos_h0);  // hour angle in radians
    D = (24.0 / pi) * H0;      // convert to hours
  }

  return D / 24.0; // fractional day length
}

phydro::PHydroResultNitrogen leaf_assimilation_rate(double fipar, double fapar,
                                                    double PAR, double PAR_max, double TAir, double VPD, double Precip, double CO2, double Nitrogen, double PA, double SWP,
                                                    double TAir_assim, double PAR_assim, double VPD_assim, double CO2_assim, double SWP_assim, double PA_assim,
                                                    phydro_canopy_parameters par, double zeta, int day){

  double infrastructure = par.infra_translation * zeta;
  phydro::ParCostNitrogen par_cost(par.alpha, par.gamma, infrastructure);
  phydro::ParPlant par_plant(par.K_leaf, par.p50_leaf, par.b_leaf);
  phydro::ParControl par_control;

  par_control.gs_method = phydro::GS_APX;
  par_control.et_method = phydro::ET_DIFFUSION;

  // day is 0-indexed within the simulated year (day 0 = Jan 1); the solar
  // formula expects 1-indexed day-of-year (N=1 for Jan 1).
  double f_day_length = day_length_fraction(HYDE_LATITUDE_DEG, day + 1);

  // TODO: this originally had the difference between acclim and inst, need to work out how to do that in CASSIA
  double Iabs_acclim = fipar * PAR_max;
  double Iabs_day    = fipar * PAR / f_day_length;
  double Iabs_24hr   = fipar * PAR;

  // TODO: check that the climate isn't altered before it is put in phydro!
  auto out_phydro_acclim = phydro::phydro_nitrogen(
    TAir_assim,     // current temperature
    TAir_assim,     // growth temperature TODO: what is this?
    Iabs_acclim,          // midday incident PAR [umol m-2 s-1]
    250,     // Net radiation [W m-2] (only used for LE calculations which we dont use) // FIXME. Should this be Rnl? See message to Beni
    VPD_assim,     // vpd [kPa]
    CO2_assim,	   // co2 [ppm]
    PA_assim,      // surface pressure [Pa]
    Nitrogen, // Leaf nitrogen store!
    fapar,                // fraction of absorbed PAR
    par.kphio,            // phi0 - quantum yield
    SWP_assim,    // soil water potential [MPa]
    par.rd,               // ratio or dark respiration to vcmax
    3,  // wind speed [m s-1], only used by PML, which we dont use, so set to global average of 3 m/s
    par.a_jmax, // TODO: a_jmax parameter
    par_plant,            // plant hydraulic traits
    par_cost,             // cost params
    par_control           // configuration params for phydro
  );

  auto photo_leaf = phydro::phydro_instantaneous_nitrogen(
    out_phydro_acclim.vcmax25, // acclimated vcmax25
    out_phydro_acclim.jmax25,  // acclimated jmax25
    TAir,            // current temperature
    TAir_assim,          // growth temperature
    Iabs_day,                  // daytime mean incident PAR [umol m-2 s-1]
    250,            // mean net radiation [W m-2] (only used for LE calculations which we dont use)
    VPD,           // vpd [kPa]
    CO2,	       // co2 [ppm]
    PA,            // surface pressure [Pa]
    fapar,                     // fraction of absorbed PAR
    par.kphio,                 // phi0 - quantum yield
    SWP,           // soil water potential [MPa]
    par.rd,                    // ratio or dark respiration to vcmax
    3,         // wind speed [m s-1], only used by PML, which we dont use, so set to global average of 3 m/s
    par.a_jmax,
    par_plant,                 // plant hydraulic traits
    par_cost,                  // cost params
    par_control                // configuration params for phydro
  );

  photo_leaf.a        *= f_day_length;
  photo_leaf.e        *= f_day_length;
  photo_leaf.gs       *= f_day_length;

  return photo_leaf;
}


// **
// ** Building the canopy
// **


double q(double z, double height, phydro_canopy_parameters par){
  if (z > height || z < 0) return 0;
  else{
    double zHn_1 = pow(z / height, par.n - 1);
    double zHn   = zHn_1 * z / height;
    return par.m * par.n * pow(1 - zHn, par.m - 1) * zHn_1;
  }
}

double zm(phydro_canopy_parameters par, double height){
  return par.zm_H * height;
}

double crown_area_above(double z, double crown_area, double height, phydro_canopy_parameters par){
  if (z == 0) return crown_area; // shortcut because z=0 is used often

  double fq = q(z, height, par) / par.qm;
  if (z >= zm(par, height)){
    return crown_area * fq * fq * (1 - par.fg);
  }
  else {
    return crown_area * (1 - fq * fq * par.fg);
  }
}

void set_forcing_acclim(double TAir, double PAR, double VPD, double CO2, double SWP, double PA,
                        double& TAir_assim, double& PAR_assim, double& VPD_assim, double& CO2_assim, double& SWP_assim, double& PA_assim,
                        phydro_canopy_parameters par){
   // Alpha, as dt is a constant of 1 day then dt = 1/365 or 1/366 as in the ecoevolutonary code
   double alpha = 1 - exp(-par.dt/par.tau_weather);
   TAir_assim += alpha * (TAir - TAir_assim);
   PAR_assim += alpha * (PAR - PAR_assim);
   VPD_assim += alpha * (VPD - VPD_assim);
   CO2_assim += alpha * (CO2 - CO2_assim);
   SWP_assim += alpha * (SWP - SWP_assim);
   PA_assim += alpha * (PA - PA_assim);
}


// **
// ** Plant assimilation for the canopy
// **


PlantAssimilationResult calc_plant_assimilation_rate(double PAR, double PAR_max, double TAir, double VPD, double Precip, double CO2, double Nitrogen, double PA, double SWP,
                                                     phydro_canopy_parameters par, double lai, double crown_area, double height, double zeta, int day){
  static double TAir_assim, PAR_assim, VPD_assim, CO2_assim, SWP_assim, PA_assim;
  if (day == 0) {
    TAir_assim = TAir;
    PAR_assim = PAR;
    VPD_assim = VPD;
    CO2_assim = CO2;
    SWP_assim = SWP;
    PA_assim = PA;
  } else {
    set_forcing_acclim(TAir, PAR, VPD, CO2, SWP, PA,
                       TAir_assim, PAR_assim, VPD_assim, CO2_assim, SWP_assim, PA_assim,
                       par);
  }

  /*
   * Calculate the parameters
   */
  // crown_area_above is the leaved area not considering gaps
  double total_crown_area = crown_area_above(0.0, crown_area, height, par); // This is just one tree rather than all of them as in PlantFate
  double n_layers = 3; // TODO: make this work int(total_crown_area / 0.99);

  // Previously hardcoded here as `lai = 1.8;` ("TODO: so this the same as
  // the amazon values") which silently discarded whatever the caller passed
  // in. Removed 2026-10-29 - now uses the real `lai` parameter (CASSIA.cpp's
  // phydro branch passes the Hyytiala-calibrated value, 4, not Amazon's 1.8).
  double fapar = 1 - exp(-par.k_light * lai);
  bool by_layer = false;

  /*
   * plant_assim processes start
   */
  PlantAssimilationResult plant_assim;
  plant_assim.gpp        = 0;
  plant_assim.rleaf      = 0;
  plant_assim.trans      = 0;
  plant_assim.dpsi_avg   = 0;
  plant_assim.vcmax_avg  = 0;
  plant_assim.vcmax25_avg = 0;
  plant_assim.mc_avg     = 0;
  plant_assim.gs_avg     = 0;
  plant_assim.c_open_avg = 0;
  plant_assim.nitrogen_avg  = 0;

  double ca_cumm = 0;
  // TODO: how do you know how many canopy layers you need? n_layers
  for (int ilayer=0; ilayer <= n_layers; ++ilayer){
    double zst = par.z_star[ilayer];

    double ca_layer = crown_area_above(zst, crown_area, height, par) - ca_cumm;

    if (by_layer == true){
      auto res = leaf_assimilation_rate(par.canopy_openness[ilayer], fapar,
                                        PAR, PAR_max, TAir, VPD, Precip, CO2, Nitrogen, PA, SWP,
                                        TAir_assim, PAR_assim, VPD_assim, CO2_assim, SWP_assim, PA_assim,
                                        par, zeta, day);
      plant_assim.gpp        += (res.a + res.vcmax * par.rd) * ca_layer;
      plant_assim.rleaf      += (res.vcmax * par.rd) * ca_layer;
      plant_assim.trans      += res.e * ca_layer;
      plant_assim.dpsi_avg   += res.dpsi * ca_layer;
      plant_assim.vcmax_avg  += res.vcmax * ca_layer;
      plant_assim.gs_avg     += res.gs * ca_layer;
      plant_assim.vcmax25_avg += res.vcmax25 * ca_layer;
      plant_assim.mc_avg     += res.mc * ca_layer;
      plant_assim.nitrogen_avg     += res.n_leaf * ca_layer;
    }

    plant_assim.c_open_avg += par.canopy_openness[ilayer] * ca_layer;
    ca_cumm += ca_layer;

  }
  if (!(fabs(ca_cumm / crown_area - 1) < 1e-6)) {
    Rcpp::Rcout << "DEBUG ca_cumm/crown_area assert about to fail: "
                << "day=" << day << " height=" << height << " crown_area=" << crown_area
                << " ca_cumm=" << ca_cumm << " ratio=" << (ca_cumm / crown_area) << "\n";
  }
  assert(fabs(ca_cumm / crown_area - 1) < 1e-6);
  double ca_total = crown_area;                   // total crown area
  plant_assim.c_open_avg /= ca_total;                // unitless
  if (by_layer == true){
    plant_assim.dpsi_avg   /= ca_total;                // MPa
    plant_assim.vcmax_avg  /= ca_total;                // umol CO2/m2/s
    plant_assim.gs_avg     /= ca_total;                // mol CO2/m2/s
    plant_assim.vcmax25_avg /= ca_total;               // umol CO2/m2/s
    plant_assim.mc_avg     /= ca_total;                // unitless
    plant_assim.nitrogen_avg     /= ca_total;                 // TODO: add units
  }

  if (by_layer == false){
    auto res = leaf_assimilation_rate(plant_assim.c_open_avg, fapar,
                                      PAR, PAR_max, TAir, VPD, Precip, CO2, Nitrogen, PA, SWP,
                                      TAir_assim, PAR_assim, VPD_assim, CO2_assim, SWP_assim, PA_assim,
                                      par, zeta, day);
    // gpp is reported to the rest of CASSIA (photosynthesis.GPP) as a per-area
    // rate in g C/m2/day, matching the PRELES/photosynthesis_as_input
    // convention - so it must NOT be multiplied by ca_total here (that would
    // make it a whole-tree absolute total, which silently scales with tree
    // size/crown growth - confirmed as the cause of GPP escalating ~3.5x
    // across a 3-year run despite near-flat weather). Confirmed against the
    // Plant-FATE/TreeLife companion codebase (src/plant.tpp divides this
    // same field by crown_area to get a per-area derivative, i.e. it is NOT
    // already per-area before that division). rleaf has no downstream
    // consumer in CASSIA (grep confirms), so its whole-tree-total scaling is
    // left as-is rather than guessing at an unused convention.
    plant_assim.gpp        = (res.a + res.vcmax * par.rd);
    plant_assim.rleaf      = (res.vcmax * par.rd) * ca_total;
    plant_assim.trans      = res.e * ca_total;
    plant_assim.dpsi_avg   = res.dpsi;
    plant_assim.vcmax_avg  = res.vcmax;
    plant_assim.gs_avg     = res.gs;
    plant_assim.vcmax25_avg = res.vcmax25;
    plant_assim.mc_avg     = res.mc;
    plant_assim.nitrogen_avg     = res.n_leaf;
  }

  // g C/m2/day, not kg biomass/unit_t: par.cbio (2.45e-2 kg biomass/mol CO2)
  // bakes in a ~2.04x carbon-to-biomass inflation factor that does not
  // belong in a quantity reported and consumed elsewhere as g C. Use the
  // direct carbon-mass conversion instead: umol CO2/m2/s -> g C/m2/day.
  plant_assim.gpp   *= (86400 * 1e-6 * 12.011);
  plant_assim.rleaf *= (86400 * 1e-6 * par.cbio);        // umol co2/s ----> umol co2/unit_t --> mol co2/unit_t --> kg/unit_t
  plant_assim.trans *= (86400 * 18e-3);                  // mol h2o/s  ----> mol h2o/unit_t  --> kg h2o /unit_t

  return plant_assim;
}

