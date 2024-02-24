#include "CASSIA.h"

struct fS_out {
  double fS;
  double S;
};

fS_out fS_model(double S, double T, p2 GPP_par) {
  double fS;

  S = S + (T-S)/GPP_par.tau;
  if (0 > S-GPP_par.S0) fS=0; else fS= S-GPP_par.S0;
  if (1 < fS/GPP_par.Smax) fS=1; else fS=fS/GPP_par.Smax;

  fS_out out;
  out.fS = fS;
  out.S = S;
  return(out);
};


/*
 * Seasonality of foliage
 */

struct fPheno_out {
  double fPheno;
  double PhenoS;
};

fPheno_out fPheno_model(p2 GPP_par, double T, double PhenoS,
                    int DOY, double fS) {
  double m;
  double fPheno=0;

  if (GPP_par.t0 > -998) { // ie not -999
    // Budbreak must occur between specified min. date and end of July
    if ( (DOY > (GPP_par.t0 - 0.5)) & (DOY < 213) )  {
      m = (T - GPP_par.tcrit);
      if (m < 0) m = 0;
      PhenoS = PhenoS + m ;
    } else {
      PhenoS = 0;
    }

    if (PhenoS > GPP_par.tsumcrit - 0.005) fPheno = 1; else fPheno = 0;
    /* Quick solution to leaf out:
     * After end of July we just apply season prediction based on conifer fS
     *  for gradual leaf out. Assume leaves drop much faster that fS.
     *  ...essentially this should be light driven process...i think. */
    if (DOY > 212) {
      fPheno = fS * fS;
      if (fPheno < 0.5) {
        fPheno = 0;
      }
    }

    /* If there is no t0 parameter, it is an evergreen */
  } else {
    fPheno = 1;
  }

  fPheno_out out;
  out.PhenoS = PhenoS;
  out.fPheno = fPheno;
  return(out);
};

struct snow_out {
  double snow;
  double snowmelt;
  double precip;
};

/*
 * Snow model
 */

snow_out snow(double TAir, double Precip, double snow, p4 SnowRain_par, double snow_melt){

  double NewSnow;
  if (TAir < SnowRain_par.SnowThreshold) {
    NewSnow=Precip;
    Precip = 0;
  } else {
    NewSnow=0;
  }

  double SnowMelt;
  if (TAir > SnowRain_par.T_0)
    SnowMelt = SnowRain_par.MeltCoef*(TAir-SnowRain_par.T_0);
  else SnowMelt=0;

  if (snow + NewSnow - SnowMelt < 0) {
    SnowMelt=NewSnow + snow;
    snow =0;
  } else {
    snow = snow + NewSnow - SnowMelt;
  }

  snow_out out;
  out.snow = snow;
  out.precip = Precip;
  out.snowmelt = SnowMelt;

  return(out);
};

/*
 * CO2
 */

double fCO2_model_mean(double CO2, p2 GPP_par ) {
  return(1 + GPP_par.bCO2 * log(CO2/380) );
};

double fCO2_ET_model_mean(double CO2, p2 GPP_par ) {
  return(1 + GPP_par.xCO2 * log(CO2/380) );
};

/*
 * Nitrogen
 */

double fN_sub(double N, p7 N_par) {
  double out = N_par.a*N+N_par.b;
  return(out);
};

/*
 * Rain interception
 */

struct interception_out {
  double precip;
  double intercepted;
};

interception_out interceptionfun(double Precip, double intercepted, double TAir,
                      p4 SnowRain_par, double fAPAR) {

  if (TAir > SnowRain_par.SnowThreshold)  {
    intercepted = Precip * (SnowRain_par.I0 * fAPAR / 0.75);
    Precip = Precip - intercepted;
  } else {
    intercepted = 0;
  }

  interception_out out;
  out.intercepted = intercepted;
  out.precip = Precip;
  return(out);
}

/*
 * ET function
 */

struct ETfun_out {
  double ET;
  double canw;
  double fE;
  double transp;
  double evap;
  double fWE;
};

ETfun_out ETfun(double D, double theta, double ppfd, double fAPAR, double T,
                p3 ET_par, p1 Site_par,
                double canw,
                double fE, double A,
                double fWgpp,  p2 GPP_par,  //double fCO2mean,
                double CO2,
                int etmodel, double transp,
                double evap, double fWE) {

  double fCO2_ET_model_mean(double CO2, p2 GPP_par);

  double thetavol = theta/Site_par.soildepth;
  double REW=(thetavol-Site_par.ThetaPWP)/(Site_par.ThetaFC-Site_par.ThetaPWP);
  //  double fEsub = -999; // Minimum of fW and fD returned if ET-model
  //			* flag indicates similar modifier as for GPP
  double fWsub=1;
  //  double fDsub=1;
  double et;
  double lambda, psychom, s; //, rho;
  double cp = 1003.5; // J/(kg K) (nearly constant, this is dry air on sea level)
  double MWratio = 0.622; // Ratio of molecular weigths of water vapor and dry air;
  // double R = 287.058; // J/(kg K) Specific gas constant for dry air, wiki
  // double zh, zm, d, zom, zoh;
  // If pressure is not inputted use default
  double pressure = 101300; // Pa

  double fCO2mean = fCO2_ET_model_mean(CO2, GPP_par);

  // rho=pressure/(R * (T+273.15) ); // Dry air density, kg/m3
  lambda = (-0.0000614342 * pow(T, 3) + 0.00158927 * pow(T, 2) -
    2.36418 * T +  2500.79) * 1000; // J/kg
  psychom= cp * pressure / (lambda* MWratio); // Pa/C, wiki
  s = 1000 * 4098.0 * (0.6109 * exp((17.27 * T)/(T+237.3))) /
    pow(T+237.3, 2);  // Pa/C! (Ice has nearly the same slope)


  // Calculate soil constraint, simple linear following Granier 1987
  if (ET_par.soilthres < -998) { /*-999 omits water control*/
  fWsub = 1;
  } else {
    if (REW < ET_par.soilthres) {
      if (REW > 0.01) fWsub = REW/ET_par.soilthres; else fWsub = 0.0;
    } else {
      fWsub = 1.0;
    }
  }


  if (canw > 0.00000001) fWsub = 1;

  fE = fWsub;
  fWE = fWsub;

  if (D < 0.01) D=0.01;

  if (etmodel == -1) {
    transp = D * ET_par.beta*A/pow(D, ET_par.kappa) *
            pow(fWgpp, ET_par.nu) * // ET differently sensitive to soil water than GPP
            fCO2mean;
    evap = ET_par.chi *  (1-fAPAR) *  fWsub * ppfd;
    et = (transp + evap) * s / (s + psychom);
  }

  if (etmodel == 0) {
    transp = D * ET_par.beta*A/pow(D, ET_par.kappa) *
            pow(fWgpp, ET_par.nu) * // ET differently sensitive to soil water than GPP
            fCO2mean;
    evap = ET_par.chi *  s / (s + psychom) * (1-fAPAR) *  fWsub * ppfd;
    et = transp + evap;
  }
  if (etmodel == 1) {
    transp = D * ET_par.beta*A/pow(D, ET_par.kappa) *
            pow(fWgpp, ET_par.nu) * // ET differently sensitive to soil water than GPP
            fCO2mean;
    evap = ET_par.chi * (1-fAPAR) *  fWsub * ppfd;
    et = transp + evap;
  }
  if (etmodel == 2) {
    et = D * (1 + ET_par.beta/pow(D, ET_par.kappa)) * A / CO2 *
        pow(fWgpp, ET_par.nu) * // ET differently sensitive to soil water than GPP
        fCO2mean +  // Mean effect of CO2 on transpiration
        ET_par.chi * (1-fAPAR) *  fWsub * ppfd;
  }

  ETfun_out out;
  out.canw = canw;
  out.ET = et;
  out.evap = evap;
  out.fE = fE;
  out.fWE = fWE;
  out.transp = transp;
  return(out);
}

/*
 * Swbalance
 */


struct sw_balance_out {
  double theta;
  double drainage;
  double snow;
  double canw;
  double theta_snow;
};

sw_balance_out swbalance(double theta, double throughfall, double snowmelt, double et,
                         p1 sitepar, double drainage,
                         double snow, double canw, p4 SnowRain_par) {

  double st0, etfromvegandsoil=0;

  // Evaporate first from wet canopy and snow on ground

  if (SnowRain_par.CWmax > 0.00000001) {
    if ( (canw + snow - et) > 0 ) {
      if ( (canw - et) > 0 ) {
        canw = canw - et;
        etfromvegandsoil = 0;
      } else if (canw - et < 0) { // in this case, there's enough snow left
        snow = snow + canw - et;
        canw = 0;
        etfromvegandsoil = 0;
      }
    } else {
      etfromvegandsoil = et - canw - snow;
      canw=0.0;
      snow = 0.0;
    }

  } else {
    if ( (snow - et) > 0 ) {
      snow = snow - et;
      etfromvegandsoil = 0;
    } else if (snow - et < 0) { // in this case, there's enough snow left
      etfromvegandsoil = et - snow;
      snow = 0;
    } else {
      snow = 0.0;
    }
  }

  et = etfromvegandsoil;

  // Water balance without drainage
  st0 = theta + throughfall + snowmelt  - et;
  if (st0 <= 0) st0 = 0.0001;

  // Calculate what is left to drainage after partial balance update above:
  if (sitepar.tauDrainage > 0) {


    // Simple time delay drainage above FC:
    if (st0 > sitepar.ThetaFC * sitepar.soildepth) {
      drainage = (st0 - sitepar.ThetaFC * sitepar.soildepth) /
        sitepar.tauDrainage;
    } else {
      drainage = 0;
    }
    theta = st0 - drainage;
  }

  sw_balance_out out;
  out.theta = theta;
  out.drainage = drainage;
  out.snow = snow;
  out.canw = canw;

  return(out);
}


/*
 * GPP function
 */


gpp_out GPPfun(double PAR, double VPD, double CO2, double theta,
            double fAPAR, double Nitrogen, double fSsub,
            p2 GPP_par, p1 Site_par, p7 N_par) {

  double thetavol = theta/Site_par.soildepth;
  double REW=(thetavol-Site_par.ThetaPWP)/(Site_par.ThetaFC-Site_par.ThetaPWP);

  // Calculate first the reference condition (ca=380 ppm) effect
  double fDsub = exp(GPP_par.kappa * VPD);
  fDsub = fDsub > 1 ? 1 : fDsub;

  double fEsub, fWsub, fLsub;
  if (GPP_par.soilthres < -998) {
    fWsub = 1.0;
  } else {
    if (REW < GPP_par.soilthres) {
      if (REW > 0.01) fWsub = REW/GPP_par.soilthres; else fWsub = 0.0;
    } else {
      fWsub = 1.0;
    }
  }

  fLsub = 1 / (GPP_par.gamma * PAR + 1);

  if (fDsub > fWsub) fEsub = fWsub; else fEsub = fDsub;
  double fW = fWsub;
  double fD = fEsub; // TODO: why is fD not used later?

  double gpp380 = GPP_par.beta * PAR * fSsub * fLsub * fEsub * fAPAR;
  double fCO2 = fCO2_model_mean(CO2, GPP_par);
  double fN = fN_sub(Nitrogen, N_par);
  double gpp = gpp380 * fCO2 * fN;

  // TODO: output
  gpp_out out;
  out.fW = fW;
  out.fD = fD;
  out.fE = fEsub; // TODO: is this fE? Or is it in the evapotranspiration?
  out.fN = fN;
  out.fCO2 = fCO2;
  out.gpp = gpp;
  out.gpp380 = gpp380;

  return(out);
}

/*
 * Preles
 */

photosynthesis_out preles(int day,
                          double PAR, double TAir, double VPD, double Precip,
                          double CO2, double fAPAR, double Nitrogen,
                          p1 Site_par,
                          p2 GPP_par,
                          p3 ET_par,
                          p4 SnowRain_par,
                          p5 Water_par,
                          p7 N_par,
                          int etmodel,
                          double theta,
                          double theta_snow,
                          double theta_canopy,
                          double Throughfall,
                          double S_state,
                          double PhenoS,
                          double Snowmelt,
                          double intercepted,
                          double Drainage,
                          double canw,
                          double fE,
                          double transp,
                          double evap,
                          double fWE,
                          double fW,
                          double gpp380)
{

  if (PAR < -900) {std::cout << "PAR is too low\n";}
  if (TAir < -900) {std::cout << "TAir is too low\n";}
  if (VPD < 0 || VPD > 6) {std::cout << "VPD is out of bounds\n";}
  // if (Precip <    0) {std::cout << "Precipitation is too low\n";}
  if (CO2 < 0) {std::cout << "CO2 is too low\n";}
  if (Nitrogen < 0) {std::cout << "Nitrogen is too low\n";}


  // TODO: should these all be 0?

  fS_out fS = fS_model(S_state, TAir, GPP_par);
  fPheno_out fPheno = fPheno_model(GPP_par, TAir, PhenoS, day, fS.fS);
  fAPAR = fAPAR * fPheno.fPheno;

  gpp_out gpp = GPPfun(PAR, VPD, CO2, theta,
                      fAPAR, Nitrogen, fS.fS,
                      GPP_par, Site_par, N_par);

  snow_out snow_values = snow(TAir, Precip, theta_snow, SnowRain_par, Snowmelt);

  Throughfall = snow_values.precip;
  interception_out interception = interceptionfun(Throughfall, intercepted,
                                                  TAir, SnowRain_par, fAPAR);

  if (SnowRain_par.CWmax <= 0.00000001) {
    Throughfall = Throughfall + interception.intercepted; // Where is through fall from?
  } else {
    if (interception.intercepted + theta_canopy > SnowRain_par.CWmax * fAPAR) {
      Throughfall = Throughfall + interception.intercepted +
        theta_canopy - SnowRain_par.CWmax  * fAPAR;
      theta_canopy = SnowRain_par.CWmax  * fAPAR;
    } else {
      theta_canopy = interception.intercepted + theta_canopy;
    }
  }
  // TOOD: theta canopy should be considered here!

  ETfun_out ET_out = ETfun(VPD, theta, PAR, fAPAR, TAir,
                           ET_par, Site_par,
                           theta_canopy,
                           fE, // Soil water constrain on evaporation
                           gpp380,
                           fW, // soil water constrain of GPP at 380 ppm
                           GPP_par, //fCO2_ET_model_mean(CO2[i], GPP_par),
                           CO2,
                           etmodel,
                           transp,
                           evap,
                           fWE);


  sw_balance_out soilwater_balance = swbalance(theta, Throughfall, snow_values.snowmelt, ET_out.fE,
                                               Site_par, Drainage,
                                               snow_values.snow, ET_out.canw, SnowRain_par);
  // TODO: throughfall should come out of the soil water balence

  double SOG = theta_snow;
  double SW = soilwater_balance.theta;
  double Canopywater = soilwater_balance.canw;

  photosynthesis_out out;
  out.GPP = gpp.gpp;
  out.ET = ET_out.ET;
  out.SoilWater = SW;
  out.fS = fS.fS;
  out.fCO2 = gpp.fCO2;
  out.fE = gpp.fE;
  out.fW = gpp.fW;
  out.fN = gpp.fN;
  out.theta_canopy = theta_canopy;
  out.Throughfall = Throughfall;
  out.theta = soilwater_balance.theta;
  out.theta_snow = soilwater_balance.theta_snow;
  out.S_state = fS.S;
  out.PhenoS = fPheno.PhenoS;
  out.Snowmelt = snow_values.snowmelt;
  out.intercepted = interception.intercepted;
  out.Drainage = soilwater_balance.drainage;
  out.canw = ET_out.canw;
  out.fE = ET_out.fE;
  out.transp = ET_out.transp;
  out.evap = ET_out.evap;
  out.fWE = ET_out.fWE;
  out.gpp380 = gpp.gpp380;

  return(out);
}

/*
 * PRELES Wrapper
 */


// [[Rcpp::export]]
Rcpp::List preles_test_cpp(int start_year, int end_year,
                           Rcpp::DataFrame weather,
                           std::vector<double> pPREL,
                           int etmodel) {

  p1 parSite = make_p1(pPREL);
  p2 parGPP = make_p2(pPREL);
  p3 parET = make_p3(pPREL);
  p4 parSnowRain = make_p4(pPREL);
  p5 parWater = make_p5(pPREL);
  p7 parN = make_p7(pPREL);

  std::vector<double> PAR = weather["PAR"];
  std::vector<double> TAir = weather["T"];
  std::vector<double> VPD = weather["VPD"];
  std::vector<double> Precip  = weather["Rain"];
  std::vector<double> CO2 = weather["CO2"];
  std::vector<double> fAPAR = weather["fAPAR"];
  std::vector<double> Nitrogen = weather["Nitrogen"];

  std::vector<int> years;
  for (int i = start_year; i <= end_year; i++) {
    years.push_back(i);
  }

  int no_days;
  photo_out_vector photosynthesis_output;
  photosynthesis_out photosynthesis_old;
  for (int year : years) {
    if (year == 2008 | year == 2012 | year == 2016 | year == 2020) {
      no_days = 366;
    } else {
      no_days = 365;
    }

    double theta, theta_snow, theta_canopy, Throughfall, S, PhenoS,
    Snowmelt, intercepted, Drainage, canw, fE, transp, evap, fWE, fW, gpp380;

    for (int day = 0; day < no_days; day++) {
      if (day == 1 & year == start_year) {
        theta = parWater.SW; // Correct
        theta_canopy = parWater.CW; // Correct
        theta_snow = parWater.SOG; // Correct
        gpp380 = 0; // Correct
        S = parWater.S; // Correct
        PhenoS = 0; // Correct
        fE = 0;
        Throughfall = 0;
        Snowmelt = 0;
        intercepted = 0;
        Drainage = 0;
        canw = 0;
        transp = 0;
        evap = 0;
        fWE = 0;
        fW = 0;
      } else {
        theta = photosynthesis_old.theta;
        theta_canopy = photosynthesis_old.theta_canopy;
        theta_snow = photosynthesis_old.theta_snow;
        Throughfall = photosynthesis_old.Throughfall;
        S = photosynthesis_old.S_state;
        PhenoS = photosynthesis_old.PhenoS;
        Snowmelt = photosynthesis_old.Snowmelt;
        intercepted = photosynthesis_old.intercepted;
        Drainage = photosynthesis_old.Drainage;
        canw = photosynthesis_old.canw;
        fE = photosynthesis_old.fE;
        transp = photosynthesis_old.transp;
        evap = photosynthesis_old.evap;
        fWE = photosynthesis_old.fWE;
        fW = photosynthesis_old.fW;
        gpp380 = photosynthesis_old.gpp380;
      }

      photosynthesis_out photosynthesis = preles(day,
                                                 PAR[day], TAir[day], VPD[day], Precip[day],
                                                 CO2[day], fAPAR[day], Nitrogen[day],
                                                 parSite,
                                                 parGPP,
                                                 parET,
                                                 parSnowRain,
                                                 parWater,
                                                 parN,
                                                 etmodel,
                                                 theta,
                                                 theta_snow,
                                                 theta_canopy,
                                                 Throughfall,
                                                 S,
                                                 PhenoS,
                                                 Snowmelt,
                                                 intercepted,
                                                 Drainage,
                                                 canw,
                                                 fE,
                                                 transp,
                                                 evap,
                                                 fWE,
                                                 fW,
                                                 gpp380);
      photosynthesis_old = photosynthesis;

      photosynthesis_output.GPP.push_back(photosynthesis.GPP);
      photosynthesis_output.ET.push_back(photosynthesis.ET);
      photosynthesis_output.SoilWater.push_back(photosynthesis.SoilWater);
      photosynthesis_output.fW.push_back(photosynthesis.fW);
      photosynthesis_output.fS.push_back(photosynthesis.fS);
      photosynthesis_output.fE.push_back(photosynthesis.fE);
      photosynthesis_output.fN.push_back(photosynthesis.fN);
    }
  }

  return Rcpp::List::create(Rcpp::_["GPP"] = photosynthesis_output.GPP,
                            Rcpp::_["ET"] = photosynthesis_output.ET,
                            Rcpp::_["SWC"] = photosynthesis_output.SoilWater,
                            Rcpp::_["fW"] = photosynthesis_output.fW,
                            Rcpp::_["fS"] = photosynthesis_output.fS,
                            Rcpp::_["fE"] = photosynthesis_output.fE,
                            Rcpp::_["fN"] = photosynthesis_output.fN);
}
