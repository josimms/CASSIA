#include "CASSIA.h"

// Vector to symphony_parameters
symphony_parameters vector_to_symphony(std::vector<double> input) {
  symphony_parameters params;

  params.A = input[0];
  params.s = input[1];
  params.r = input[2];
  params.m_p = input[3];
  params.e_p = input[4];
  params.r_p = input[5];
  params.k = input[6];
  params.alpha = input[7];
  params.beta = input[8];
  params.i = input[9];
  params.l = input[10];
  params.y = input[11];
  params.u = input[12];
  params.e = input[13];
  params.Ca = input[14];
  params.psi_i = input[15];

  return(params);
}

// [[Rcpp::export]]
Rcpp::List symphony(std::vector<double> params) {
  // Parameters into symphony_parameters format
  symphony_parameters symphony_params = vector_to_symphony(params);

  // Make vectors
  std::vector<double> psi_ph;
  std::vector<double> psi_up;
  std::vector<double> psi_ims;
  std::vector<double> psi_imf;
  std::vector<double> psi_d;
  std::vector<double> psi_f;
  std::vector<double> C_p;
  std::vector<double> C_ds;
  std::vector<double> C_df;
  std::vector<double> C_f;
  std::vector<double> N;
  std::vector<double> C_s;

  double psi_ph_min;
  double psi_d_min;
  double psi_f_min;

  // Initialisation with the observed values of the paper
  // TODO: change the values, so they are not har coded and so they are boreal
  C_p.push_back(525);
  C_f.push_back(635);
  N.push_back(2.09);
  // Initialised with 1 to test the code
  psi_ph.push_back(1);
  psi_up.push_back(1);
  psi_ims.push_back(1);
  psi_imf.push_back(1);
  psi_f.push_back(1);
  psi_d.push_back(1);

  C_ds.push_back(1);
  C_df.push_back(1);
  C_s.push_back(1);

  // TODO: add the initial conditions and decide which process should happen first

  for (int day = 0; day <= 364; day++) { // Indexes start from the 2nd day of the year
    // Photosynthesis
    psi_ph_min = std::min(symphony_params.k * symphony_params.Ca, (symphony_params.e * N[day])/symphony_params.beta + symphony_params.r_p * C_p[day]);
    psi_ph.push_back(psi_ph_min);
    // N in
    psi_up.push_back(symphony_params.beta * (psi_ph[day] - symphony_params.r_p*C_p[day]));
    // N immobilisation, mineralisation
    psi_ims.push_back(symphony_params.alpha * symphony_params.r * C_ds[day] + (symphony_params.beta - symphony_params.alpha)*psi_d[day]);
    // N immobilisation,
    psi_imf.push_back(symphony_params.alpha * symphony_params.r * C_df[day] + (symphony_params.beta - symphony_params.alpha)*psi_f[day]);
    // N decomposition flux
    psi_d_min = std::min(symphony_params.y*C_f[day], (symphony_params.i*N[day] + symphony_params.alpha * symphony_params.r * C_ds[day])/(symphony_params.alpha - symphony_params.beta));
    psi_d.push_back(psi_d_min);
    psi_f_min = std::min(symphony_params.u*C_f[day], (symphony_params.i*N[day] + symphony_params.alpha * symphony_params.r * C_df[day])/(symphony_params.alpha - symphony_params.beta));
    psi_f.push_back(psi_d_min);

    C_p.push_back(C_p[day] + psi_ph[day] - C_p[day]*(symphony_params.r_p + symphony_params.m_p + symphony_params.e_p));
    C_f.push_back(C_f[day] + symphony_params.m_p * C_p[day] - psi_d[day] - psi_f[day]);
    C_df.push_back(C_df[day] + psi_f[day] - (symphony_params.s + symphony_params.r) * C_df[day]);
    C_ds.push_back(C_ds[day] + (symphony_params.A - symphony_params.s - symphony_params.r)*C_ds[day] + psi_d[day]);
    N.push_back(N[day] + symphony_params.psi_i - symphony_params.l*N[day] - psi_up[day] + psi_ims[day] + psi_imf[day]);
    C_s.push_back(C_s[day] + (symphony_params.s - symphony_params.A) * C_ds[day]);
  }

  // Output
  return Rcpp::List::create(Rcpp::_["C_df"] = C_df,
                            Rcpp::_["C_ds"] = C_ds,
                            Rcpp::_["C_p"] = C_p,
                            Rcpp::_["C_f"] = C_f,
                            Rcpp::_["N"] = N,
                            Rcpp::_["C_s"] = C_s);
}


// [[Rcpp::export]]
Rcpp::List symphony_plus(std::vector<double> params, std::vector<double> Photosynthesis) {
  // Parameters into symphony_parameters format
  symphony_parameters symphony_params = vector_to_symphony(params);

  // Make vectors
  std::vector<double> psi_ph;
  std::vector<double> psi_up;
  std::vector<double> psi_ims;
  std::vector<double> psi_imf;
  std::vector<double> psi_d;
  std::vector<double> psi_f;
  std::vector<double> C_p;
  std::vector<double> C_ds;
  std::vector<double> C_df;
  std::vector<double> C_f;
  std::vector<double> N;
  std::vector<double> C_s;

  double psi_ph_min;
  double psi_d_min;
  double psi_f_min;

  // Initialisation with the observed values of the paper
  C_p.push_back(525);
  C_f.push_back(635);
  N.push_back(2.09);
  // PRELES ultermate input
  psi_ph = Photosynthesis;
  // Initialized with 1 to test the code
  psi_up.push_back(1);
  psi_ims.push_back(1);
  psi_imf.push_back(1);
  psi_f.push_back(1);
  psi_d.push_back(1);
  C_ds.push_back(1);
  C_df.push_back(1);
  C_s.push_back(1);

  // alpha and beta should maybe be reset each year, so I need to add them!

  // TODO: add the initial conditions and decide which process should happen first

  for (int day = 0; day <= 364; day++) { // Indexes start from the 2nd day of the year
    // N in
    psi_up.push_back(symphony_params.beta * (psi_ph[day] - symphony_params.r_p*C_p[day]));
    // N immobilisation, mineralisation
    psi_ims.push_back(symphony_params.alpha * symphony_params.r * C_ds[day] + (symphony_params.beta - symphony_params.alpha)*psi_d[day]);
    // N immobilisation,
    psi_imf.push_back(symphony_params.alpha * symphony_params.r * C_df[day] + (symphony_params.beta - symphony_params.alpha)*psi_f[day]);
    // N decomposition flux
    psi_d_min = std::min(symphony_params.y*C_f[day], (symphony_params.i*N[day] + symphony_params.alpha * symphony_params.r * C_ds[day])/(symphony_params.alpha - symphony_params.beta));
    psi_d.push_back(psi_d_min);
    psi_f_min = std::min(symphony_params.u*C_f[day], (symphony_params.i*N[day] + symphony_params.alpha * symphony_params.r * C_df[day])/(symphony_params.alpha - symphony_params.beta));
    psi_f.push_back(psi_d_min);

    C_p.push_back(C_p[day] + psi_ph[day] - C_p[day]*(symphony_params.r_p + symphony_params.m_p + symphony_params.e_p));
    C_f.push_back(C_f[day] + symphony_params.m_p * C_p[day] - psi_d[day] - psi_f[day]);
    C_df.push_back(C_df[day] + psi_f[day] - (symphony_params.s + symphony_params.r) * C_df[day]);
    C_ds.push_back(C_ds[day] + (symphony_params.A - symphony_params.s - symphony_params.r)*C_ds[day] + psi_d[day]);
    N.push_back(N[day] + symphony_params.psi_i - symphony_params.l*N[day] - psi_up[day] + psi_ims[day] + psi_imf[day]);
    C_s.push_back(C_s[day] + (symphony_params.s - symphony_params.A) * C_ds[day]);
  }

  // Output
  return Rcpp::List::create(Rcpp::_["C_df"] = C_df,
                            Rcpp::_["C_ds"] = C_ds,
                            Rcpp::_["C_p"] = C_p,
                            Rcpp::_["C_f"] = C_f,
                            Rcpp::_["N"] = N,
                            Rcpp::_["C_s"] = C_s);
}

// [[Rcpp::export]]
Rcpp::List symphony_plus_daily(std::vector<double> params,
                               double Photosynthesis, // TODO: check the units here!
                               double C_plant, // TODO: this should come from CASSIA
                               double C_FOM, // TODO: this will change
                               double N) // TODO: should be from different types of N)
{
  // Parameters into symphony_parameters format
  symphony_parameters symphony_params = vector_to_symphony(params);

  // Make values
  // Initialized with 1 to test the code
  double psi_up = 1; // TODO: work out what this is
  double psi_ims = 1;
  double psi_imf = 1;
  double psi_f = 1;
  double psi_d = 1;
  double C_decompose_SOM = 1;
  double C_decompose_FOM = 1;
  double C_SOM = 1;

  double Photosynthesis_min;
  double psi_d_min;
  double psi_f_min;

  // alpha and beta should maybe be reset each year, so I need to add them!

  // TODO: add the initial conditions and decide which process should happen first

  // N in
  psi_up = symphony_params.beta * (Photosynthesis - symphony_params.r_p*C_plant);
  // N immobilisation, mineralisation
  psi_ims = symphony_params.alpha * symphony_params.r * C_decompose_SOM + (symphony_params.beta - symphony_params.alpha)*psi_d;
  // N immobilisation,
  psi_imf = symphony_params.alpha * symphony_params.r * C_decompose_FOM + (symphony_params.beta - symphony_params.alpha)*psi_f;
  // N decomposition flux
  psi_d_min = std::min(symphony_params.y*C_FOM, (symphony_params.i*N + symphony_params.alpha * symphony_params.r * C_decompose_SOM)/(symphony_params.alpha - symphony_params.beta));
  psi_d = psi_d_min;
  psi_f_min = std::min(symphony_params.u*C_FOM, (symphony_params.i*N + symphony_params.alpha * symphony_params.r * C_decompose_FOM)/(symphony_params.alpha - symphony_params.beta));
  psi_f = psi_d_min;

  C_plant = C_plant + Photosynthesis - C_plant*(symphony_params.r_p + symphony_params.m_p + symphony_params.e_p);
  C_FOM = C_FOM + symphony_params.m_p * C_plant - psi_d - psi_f;
  C_decompose_FOM = C_decompose_FOM + psi_f - (symphony_params.s + symphony_params.r) * C_decompose_FOM;
  C_decompose_SOM = C_decompose_SOM + (symphony_params.A - symphony_params.s - symphony_params.r)*C_decompose_SOM + psi_d;
  N = N + symphony_params.psi_i - symphony_params.l*N - psi_up + psi_ims + psi_imf;
  C_SOM = C_SOM + (symphony_params.s - symphony_params.A) * C_decompose_SOM;

  // Output
  return Rcpp::List::create(Rcpp::_["C_decompose_FOM"] = C_decompose_FOM,
                            Rcpp::_["C_decompose_SOM"] = C_decompose_SOM,
                            Rcpp::_["C_plant"] = C_plant,
                            Rcpp::_["C_FOM"] = C_FOM,
                            Rcpp::_["N"] = N,
                            Rcpp::_["C_SOM"] = C_SOM);
}


