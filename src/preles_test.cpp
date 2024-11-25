#include "CASSIA.h"

photosynthesis_out preles_cpp(int day,
                              double I,
                              double T,
                              double P,
                              double D,
                              double CO2,
                              double fAPAR,
                              p1 Site_par,
                              p2 GPP_par,
                              p3 ET_par,
                              p4 SnowRain_par,
                              p5 Initials_snow,
                              double LOGFLAG) {
  // TODO: initialise / make these into values that are inputs as well!
  static double S_state, PhenoS, fD, fW, fEgpp, gpp380;
  static double theta, theta_canopy, theta_snow, Snowmelt;
  static double Interception, Drainage, etmodel, transp, evap;
  static double GPP, fE, fWE, fPheno;

  photosynthesis_out out;
  if (day==0) {
    theta=Initials_snow.SW;
    theta_canopy=Initials_snow.CW;
    theta_snow=Initials_snow.SOG;
    S_state = Initials_snow.S;
    Snowmelt = 0.0;
    Interception = 0.0;
    Drainage = 0.0;
    transp = 0.0;
    evap = 0.0;
    fWE = 0.0;
    PhenoS=0.0;
    fPheno=0.0;
    GPP = 0.0;
    gpp380 = 0.0;
    fD = 0.0;
    fW = 0.0;
    fEgpp = 0.0;
  }

  double fS = fS_model(&S_state, T, GPP_par);

  fPheno = fPheno_model(GPP_par, T, &PhenoS, day, fS);

  fAPAR = fAPAR * fPheno;

  GPPfun(&GPP, &gpp380, I, D, CO2, theta, fAPAR, fS,
         GPP_par, Site_par,  &fD, &fW, &fEgpp,
         LOGFLAG );

  Snow(T, &P, &theta_snow, SnowRain_par, &Snowmelt);

  double Throughfall = P;
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

  double ET = ETfun(D, theta, I, fAPAR, T,
                    ET_par, Site_par,
                    &theta_canopy,
                    &fE, // Soil water constrain on evaporation
                    gpp380,
                    fW, // soil water constrain of GPP at 380 ppm
                    GPP_par, //fCO2_ET_model_mean(CO2[i], GPP_par),
                    CO2,
                    LOGFLAG, etmodel,
                    &transp,
                    &evap, &fWE);

  swbalance(&theta, Throughfall, Snowmelt, ET,
            Site_par, &Drainage,
            &theta_snow, &theta_canopy, SnowRain_par);

  out.GPP = GPP;
  out.ET = ET;
  out.SoilWater = theta;
  out.fS = fS;

  return out;
}

// [[Rcpp::export]]
Rcpp::DataFrame preles_test(Rcpp::DataFrame weather) {
  std::vector<double> PAR = weather["PAR"];
  std::vector<double> VPD = weather["VPD"];
  std::vector<double> fAPAR = weather["fAPAR"];
  std::vector<double> TAir = weather["T"];
  std::vector<double> Precip  = weather["Rain"];
  std::vector<double> CO2 = weather["CO2"];

  p1 Site_par;
  Site_par.soildepth = 413.0;
  Site_par.ThetaFC = 0.45;
  Site_par.ThetaPWP = 0.118;
  Site_par.tauDrainage = 3;

  p2 GPP_par;
  GPP_par.beta = 0.745700;
  GPP_par.tau = 10.930000;
  GPP_par.S0 = -3.06300;
  GPP_par.Smax = 17.720000;
  GPP_par.kappa = -0.102700;
  GPP_par.gamma = 0.036730;
  GPP_par.soilthres = 0.777900;
  GPP_par.bCO2 = 0.5;
  GPP_par.xCO2 = -0.364000;
  GPP_par.t0 = -999;
  GPP_par.tcrit = -999;
  GPP_par.tsumcrit = -999;

  p3 ET_par;
  ET_par.beta = 0.271500;
  ET_par.kappa = 0.835100;
  ET_par.chi = 0.073480;
  ET_par.soilthres = 0.999600;
  ET_par.nu = 0.442800;

  p4 SnowRain_par;
  SnowRain_par.MeltCoef = 1.2;
  SnowRain_par.I0 = 0.33;
  SnowRain_par.CWmax = 4.970496;
  SnowRain_par.SnowThreshold = 0.0;
  SnowRain_par.T_0 = 0.0;

  p5 Initials_snow;
  Initials_snow.SW = 160,0;
  Initials_snow.CW = 0.0;
  Initials_snow.SOG = 0.0;
  Initials_snow.S = 20.0;

  double LOGFLAG = 0;

  int days_per_year = PAR.size();
  photosynthesis_out out;
  photo_out_vector out_vector;

  for (int day = 0; day < days_per_year; day++) {

    double I = PAR[day];
    double T = TAir[day];
    double P = Precip[day];
    double D = VPD[day];
    double co2 = CO2[day];
    double fapar = fAPAR[day];

    out = preles_cpp(day, I, T, P, D, co2, fapar,
                     Site_par, GPP_par, ET_par, SnowRain_par,
                     Initials_snow, LOGFLAG);

    out_vector.GPP.push_back(out.GPP);
    out_vector.ET.push_back(out.ET);
    out_vector.SoilWater.push_back(out.SoilWater);
  }

  Rcpp::DataFrame out_r =  Rcpp::DataFrame::create(
    Rcpp::Named("GPP") = out_vector.GPP,
    Rcpp::Named("ET") = out_vector.ET,
    Rcpp::Named("SoilWater") = out_vector.SoilWater
  );

  return out_r;
}
