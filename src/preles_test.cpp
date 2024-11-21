#include "CASSIA.h"

/*
 * Functions not defined in the header file so should be defined here!
 */

// Declare the external C function
extern "C" {
  double fS_model(double* S, double T, p2 GPP_par);
}


/*
 * PRELES Wrapper from C++
 */

double preles_day(const std::vector<double>& S, double T, const p2& GPP_par) {
  /// Ensure the input array is not empty
  if (S.empty()) {
    throw std::invalid_argument("Input vector S is empty");
  }

  // Call the C function using the raw pointer to the vector's data
  return fS_model(const_cast<double*>(S.data()), T, GPP_par);
}


/*

fPheno = fPheno_model(GPP_par, T, &PhenoS, day, fS);

fAPAR = fAPAR * fPheno;

GPPfun(&GPP, &gpp380, I, D, CO2, theta, fAPAR, fS,
       GPP_par, Site_par,  &fD, &fW, &fEgpp,
       LOGFLAG );

// if (LOGFLAG > 1.5)
// fprintf(flog,
// "   preles(): estimated GPP=%lf\tfD=%lf\tfEgpp=%lf\n GPP380ppm %lf\n",
// GPP[i], fD[i], fEgpp, gpp380);


Snow(T, &P, &theta_snow, SnowRain_par, &Snowmelt);

// NOTE: interception model could be better
Throughfall = P;
interceptionfun(&Throughfall, &Interception, T,
                SnowRain_par, fAPAR);


if (SnowRain_par.CWmax <= 0.00000001) {
  Throughfall = Throughfall + Interception;
} else {
  if (Interception + theta_canopy > SnowRain_par.CWmax * fAPAR) {
    Throughfall = Throughfall + Interception +
      theta_canopy - SnowRain_par.CWmax  * fAPAR;
    theta_canopy = SnowRain_par.CWmax  * fAPAR;
  } else {
    theta_canopy = Interception + theta_canopy;
  }
}

ET = ETfun(D, theta, I, fAPAR, T,
           ET_par, Site_par,
           &theta_canopy,
           &fE, // Soil water constrain on evaporation
           gpp380,
           fW, // soil water constrain of GPP at 380 ppm
           GPP_par, //fCO2_ET_model_mean(CO2[i], GPP_par),
           CO2,
           etmodel,
           &transp,
           &evap, &fWE);


           swbalance(&theta, Throughfall, Snowmelt, ET,
                     Site_par, &Drainage, //&Psi[i], &Ks[i],
                     &theta_snow, &theta_canopy, SnowRain_par);

           SOG = theta_snow;
           SW = theta;
           Canopywater = theta_canopy;
           S=S_state;


           photosynthesis_out out;
           out.GPP = *GPP;
           out.ET = *ET;
           out.SoilWater = *SW;
           out.fS = *fS;
           out.fE = *fE;
           out.fW = *fW;
           // out.fN = *fN;
           out.theta_canopy = theta_canopy;
           out.Throughfall = *Throughfall;
           out.theta = theta;
           out.theta_snow = theta_snow;
           out.PhenoS = PhenoS;
           out.Snowmelt = *Snowmelt;
           out.intercepted = *Interception;
           out.Drainage = *Drainage;
           out.canw = Canopywater;
           out.transp = *transp;
           out.evap = *evap;
           out.fWE = *fWE;
           out.gpp380 = gpp380;
           return(out);
 */
