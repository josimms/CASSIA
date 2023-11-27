#include "CASSIA.h"

double fS_model(double S, double T, p2 GPP_par) {
  double fS;

  S = S + (T-S)/GPP_par.tau;
  if (0 > S-GPP_par.S0) fS=0; else fS= S-GPP_par.S0;
  if (1 < fS/GPP_par.Smax) fS=1; else fS=fS/GPP_par.Smax;

  return(fS);
};

/*
 * Seasonality of foliage
 */

double fPheno_model(p2 GPP_par, double T, double PhenoS,
                    int DOY, double fS) {
  double m;
  double fPheno=0;

  if (GPP_par.t0 > -998) { // ie not -999
    /* Budbreak must occur between specified min. date and end of July */
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

  return(fPheno);
};

struct snow_out {
  double snow;
  double snowmelt;
  double precip;
};

/*
 * Snow model
 */

snow_out snow(double TAir, double Precip, double snow, p4 SnowRain_par){

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
  out.snowmelt = SnowMelt;
  out.precip = Precip;

  return(out);
};

double fCO2_model_mean(double CO2, p2 GPP_par ) {
  return(1 + GPP_par.bCO2 * log(CO2/380) );
};

double fCO2_ET_model_mean(double CO2, p2 GPP_par ) {
  return(1 + GPP_par.xCO2 * log(CO2/380) );
};

double fN_sub(double N, p7 N_par) {
  double out = N_par.a*N+N_par.b;
  return(out);
};

/*
 * Rain interception
 */

double interceptionfun(double Precip, double TAir,
                      p4 SnowRain_par, double fAPAR) {

  double intercepted;
  if (TAir > SnowRain_par.SnowThreshold)  {
    intercepted = Precip * (SnowRain_par.I0 * fAPAR / 0.75);
    Precip = Precip - intercepted;
  } else {
    intercepted = 0;
  }

  return(intercepted);
}

/*
 * ET function
 */

double ETfun(double D, double theta, double ppfd, double fAPAR, double T,
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
  //  double fEsub = -999; /* Minimum of fW and fD returned if ET-model
  //			* flag indicates similar modifier as for GPP */
  double fWsub=1;
  //  double fDsub=1;
  double et;
  double lambda, psychom, s; //, rho;
  double cp = 1003.5; // J/(kg K) (nearly constant, this is dry air on sea level)
  double MWratio = 0.622; // Ratio of molecular weigths of water vapor and dry air;
  // double R = 287.058; // J/(kg K) Specific gas constant for dry air, wiki
  // double zh, zm, d, zom, zoh;
  /*If pressure is not inputted use default */
  double pressure = 101300; // Pa

  double fCO2mean = fCO2_ET_model_mean(CO2, GPP_par);

  // rho=pressure/(R * (T+273.15) ); // Dry air density, kg/m3
  lambda = (-0.0000614342 * pow(T, 3) + 0.00158927 * pow(T, 2) -
    2.36418 * T +  2500.79) * 1000; // J/kg
  psychom= cp * pressure / (lambda* MWratio); // Pa/C, wiki
  s = 1000 * 4098.0 * (0.6109 * exp((17.27 * T)/(T+237.3))) /
    pow(T+237.3, 2);  // Pa/C! (Ice has nearly the same slope)


  /* Calculate soil constraint, simple linear following Granier 1987*/
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

  return(et);
}

/*
 * Swbalance
 */

struct sw_balance_out {
  double theta;
  double drainage;
};

sw_balance_out swbalance(double theta, double throughfall, double snowmelt, double et,
                         p1 sitepar, double drainage,
                         double snow, double canw, p4 SnowRain_par) {

  double st0, etfromvegandsoil=0;

  /* Evaporate first from wet canopy and snow on ground */

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

  /* Water balance without drainage */
  st0 = theta + throughfall + snowmelt  - et;
  if (st0 <= 0) st0 = 0.0001;

  /* Calculate what is left to drainage after partial balance update above: */
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

  return(out);
}



/*
 * GPP function
 */

double GPPfun(double PAR, double VPD, double CO2, double theta,
            double fAPAR, double Nitrogen, double fSsub,
            p2 GPP_par, p1 Site_par, p7 N_par) {

  double thetavol = theta/Site_par.soildepth;
  double REW=(thetavol-Site_par.ThetaPWP)/(Site_par.ThetaFC-Site_par.ThetaPWP);

  /* Calculate first the reference condition (ca=380 ppm) effect */
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
  double fD = fEsub;

  double gpp380 = GPP_par.beta * PAR * fSsub * fLsub * fEsub * fAPAR;
  double fCO2 = fCO2_model_mean(CO2, GPP_par);
  double fN = fN_sub(Nitrogen, N_par);
  double gpp = gpp380 * fCO2 * fN;

  return(gpp);
}

/*
 * Preles
 */

photosynthesis_out preles(int NofDays, int day,
                          double PAR, double TAir, double VPD, double Precip,
                          double CO2,
                          double fAPAR,
                          double Nitrogen,
                          p1 Site_par,
                          p2 GPP_par,
                          p3 ET_par,
                          p4 SnowRain_par,
                          p5 Water_par,
                          p7 N_par,
                          int etmodel)
{

  if (PAR < -900) {std::cout << "PAR is too low\n";}
  if (TAir < -900) {std::cout << "TAir is too low\n";}
  if (VPD < 0 || VPD > 6) {std::cout << "VPD is out of bounds\n";}
  if (Precip <    0) {std::cout << "Precipitation is too low\n";}
  if (CO2 < 0) {std::cout << "CO2 is too low\n";}
  if (Nitrogen < 0) {std::cout << "Nitrogen is too low\n";}

  double theta, theta_snow, theta_canopy, S_state;
  double PhenoS = 0;
  double fPheno = 0;
  double fEgpp = 0;
  double gpp380 = 0;
  double fE = 0;
  double fW = 0;
  double transp = 0;
  double evap = 0;
  double fWE = 0;
  double Drainage = 0;

  double fS = fS_model(S_state, TAir, GPP_par);
  fPheno = fPheno_model(GPP_par, TAir, PhenoS, day, fS);
  fAPAR = fAPAR * fPheno;

  double gpp = GPPfun(PAR, VPD, CO2, theta,
                      fAPAR, Nitrogen, fS,
                      GPP_par, Site_par, N_par);

  snow_out snow_values = snow(TAir, Precip, theta_snow, SnowRain_par);

  double Throughfall = snow_values.precip;
  double interception = interceptionfun(Throughfall, TAir, SnowRain_par, fAPAR);

  if (SnowRain_par.CWmax <= 0.00000001) {
    Throughfall = Throughfall + interception; // Where is through fall from?
  } else {
    if (interception + theta_canopy > SnowRain_par.CWmax * fAPAR) {
      Throughfall = Throughfall + interception +
        theta_canopy - SnowRain_par.CWmax  * fAPAR;
      theta_canopy = SnowRain_par.CWmax  * fAPAR;
    } else {
      theta_canopy = interception + theta_canopy;
    }
  }

  double ET = ETfun(VPD, theta, PAR, fAPAR, TAir,
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
                    fWE); // TOOD: make sure that these make sense!


  sw_balance_out soilwater_balance = swbalance(theta, Throughfall, snow_values.snowmelt, ET,
                                               Site_par, Drainage,
                                               theta_snow, theta_canopy, SnowRain_par);
  double SoilWater; // TODO: work out the output here!

  double SOG = theta_snow;
  double SW = soilwater_balance.theta;
  double Canopywater = theta_canopy;

  photosynthesis_out out;
  out.GPP = gpp;
  out.ET = ET;
  out.SoilWater = SoilWater;
  out.S = fS;

  return(out);
}


/*
 * PRELES Wrapper
 */

//[[Rcpp::export]]
Rcpp::List preles_test_cpp(int NofDays, int day,
                           Rcpp::DataFrame weather,
                           std::vector<double> pPREL,
                           int etmodel) {

  p1 parSite = make_p1(pPREL);
  p2 parGPP = make_p2(pPREL);
  p3 parET = make_p3(pPREL);
  p4 parSnowRain = make_p4(pPREL);
  p5 parWater = make_p5(pPREL);
  p7 parN = make_p7(pPREL);

  double PAR = weather["PAR"];
  double TAir = weather["T"];
  double VPD = weather["VPD"];
  double Precip  = weather["Rain"];
  double CO2 = weather["CO2"];
  double fAPAR = weather["fAPAR"];
  double Nitrogen = weather["Nitrogen"];

  photosynthesis_out out = preles(NofDays, day,
                                  PAR, TAir, VPD, Precip,
                                  CO2,
                                  fAPAR,
                                  Nitrogen,
                                  parSite,
                                  parGPP,
                                  parET,
                                  parSnowRain,
                                  parWater,
                                  parN,
                                  etmodel);

  return Rcpp::List::create(Rcpp::_["GPP"] = out.GPP,
                            Rcpp::_["ET"] = out.ET,
                            Rcpp::_["SWC"] = out.SoilWater,
                            Rcpp::_["S"] = out.S);

}
