#include "external/phydro/inst/include/phydro.h"
#include "function_structures.h"
#include "final_parameters.h"

// **
// ** Gross and Net Assimilation
// **
phydro::PHydroResultNitrogen leaf_assimilation_rate(double fipar, double fapar,
                                                    double PAR, double TAir, double VPD, double Precip, double CO2, double Nitrogen, double PA, double SWP,
                                                    phydro_canopy_parameters par, double zeta){

  double infrastructure = par.infra_translation * zeta;
  phydro::ParCostNitrogen par_cost(par.alpha, par.gamma, infrastructure);
  phydro::ParPlant par_plant(par.K_leaf, par.p50_leaf, par.b_leaf);
  phydro::ParControl par_control;

  par_control.gs_method = phydro::GS_APX;
  par_control.et_method = phydro::ET_DIFFUSION;

  double f_day_length = 0.5;

  // TODO: this originally had the difference between acclim and inst, need to work out how to do that in CASSIA
  double Iabs_acclim = fipar * PAR;
  double Iabs_day    = fipar * PAR / f_day_length;
  double Iabs_24hr   = fipar * PAR;

  // TODO: check that the climate isn't altered before it is put in phydro!
  auto out_phydro_acclim = phydro::phydro_nitrogen(
    TAir,     // current temperature
    TAir,     // growth temperature TODO: what is this?
    Iabs_acclim,          // midday incident PAR [umol m-2 s-1]
    PAR,     // Net radiation [W m-2] (only used for LE calculations which we dont use) // FIXME. Should this be Rnl? See message to Beni
    VPD,    // vpd [kPa]
    CO2,	  // co2 [ppm]
    PA,     // surface pressure [Pa]
    Nitrogen, // Leaf nitrogen store!
    fapar,                // fraction of absorbed PAR
    par.kphio,            // phi0 - quantum yield
    SWP,    // soil water potential [MPa]
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
    TAir,          // growth temperature
    Iabs_day,                  // daytime mean incident PAR [umol m-2 s-1]
    PAR,            // mean net radiation [W m-2] (only used for LE calculations which we dont use)
    VPD,           // vpd [kPa]
    CO2,	       // co2 [ppm]
    PA,            // surface pressure [Pa]
    Nitrogen, // TODO: this should be the nitrogen_store when there is one
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

// **
// ** Plant assimilation for the canopy
// **

PlantAssimilationResult calc_plant_assimilation_rate(double fipar,
                                                     double PAR, double TAir, double VPD, double Precip, double CO2, double Nitrogen, double PA, double SWP,
                                                     phydro_canopy_parameters par, double lai, double n_layers, double crown_area, double height, double zeta){
  // TODO: lai should be calculated in the CASSIA model
  double fapar = 1 - exp(-par.k_light * lai);
  bool by_layer = false;

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

    // TODO: environment - what is this?

    if (by_layer == true){
      auto res = leaf_assimilation_rate(par.canopy_openness[ilayer], fapar,
                                        PAR, TAir, VPD, Precip, CO2, Nitrogen, PA, SWP,
                                        par, zeta);
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
                                      PAR, TAir, VPD, Precip, CO2, Nitrogen, PA, SWP,
                                      par, zeta);
    plant_assim.gpp        = (res.a + res.vcmax * par.rd) * ca_total;
    plant_assim.rleaf      = (res.vcmax * par.rd) * ca_total;
    plant_assim.trans      = res.e * ca_total;
    plant_assim.dpsi_avg   = res.dpsi;
    plant_assim.vcmax_avg  = res.vcmax;
    plant_assim.gs_avg     = res.gs;
    plant_assim.vcmax25_avg = res.vcmax25;
    plant_assim.mc_avg     = res.mc;
    plant_assim.nitrogen_avg     = res.n_leaf;
  }

  plant_assim.gpp   *= (1e-6 * par.cbio);        // umol co2/s ----> umol co2/unit_t --> mol co2/unit_t --> kg/unit_t
  plant_assim.rleaf *= (1e-6 * par.cbio);        // umol co2/s ----> umol co2/unit_t --> mol co2/unit_t --> kg/unit_t
  plant_assim.trans *= (18e-3);                  // mol h2o/s  ----> mol h2o/unit_t  --> kg h2o /unit_t

  return plant_assim;
}
