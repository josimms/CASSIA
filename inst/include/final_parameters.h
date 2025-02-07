#ifndef FINAL_PARAMETERS_H
#define FINAL_PARAMETERS_H

/*
 * CASSIA Structures
 */

struct CASSIA_parameters {
  double Q10_N;
  double Rm_N;
  double Q10_S;
  double Rm_S;
  double Q10_R;
  double Rm_R;
  double sR0;
  double sRc;
  double growth_myco;
  double root_lifetime;
  double HH0;
  double sH0;
  double LH;
  double LH0;
  double sHc;
  double sN0;
  double LN;
  double LN0;
  double sNc;
  double HN0;
  double sD0_Trad;
  double LD;
  double LD0;
  double sDc;
  double sDc_T_count;
  double tau_Ee;
  double tau_El;
  double tau_We;
  double tau_Wl;
  double tau_GPP;
  double Uggla;
  double sB0;
  double sBc;
  double LB;
  double cell_d_ew;
  double cell_d_lw;
  double cell_l_ew;
  double cell_l_lw;
  double cell_wall_density_ew;
  double cell_wall_density_lw;
  double wall_thickness_ew;
  double wall_thickness_lw;
  double cell_volume_growth_per_day_ew;
  double cell_volume_growth_per_day_lw;
  double density_tree;
  double carbon_share;
  double D0;
  double h0;
  double n_age;
  double n_length;
  double h_increment;
  double SLA;
  double LR0;
  double mycorrhiza_threshold;
  double starch0;
  double sugar0;
  double starch_needles0;
  double starch_phloem0;
  double starch_xylem_sh0;
  double starch_xylem_st0;
  double starch_roots0;
  double sugar_needles0;
  double sugar_phloem0;
  double sugar_roots0;
  double sugar_xylem_sh0;
  double sugar_xylem_st0;
  double Wala_needles;
  double Wala_phloem;
  double Wala_xylem_sh;
  double Wala_xylem_st;
  double Wala_roots;
  double carbon_sugar;
  double carbon_starch;
  double alfa_needles;
  double alfa_phloem;
  double alfa_xylem_sh;
  double alfa_xylem_st;
  double alfa_roots;
  double tau_s;
  double tau_t;
  double starch00;
  double sugar00;
  double Wala;
  double alfa;
  double Q10s;
  double Q10d;
  double SCb;
  double sugar_level;
  double Ad0_needles;
  double Ad0_phloem;
  double Ad0_roots;
  double Ad0_xylem_sh;
  double Ad0_xylem_st;
  double lambda_needles;
  double lambda_phloem;
  double lambda_roots;
  double lambda_xylem_sh;
  double lambda_xylem_st;
  double delta_needles;
  double delta_phloem;
  double delta_roots;
  double delta_xylem_sh;
  double delta_xylem_st;
  double lower_bound_needles;
  double lower_bound_phloem;
  double lower_bound_roots;
  double lower_bound_xylem_sh;
  double lower_bound_xylem_st;
  double tau_emergancy_needles;
  double tau_emergancy_phloem;
  double tau_emergancy_roots;
  double tau_emergancy_xylem_sh;
  double tau_emergancy_xylem_st;
  double resistance_needles_to_phloem;
  double resistance_phloem_to_roots;
  double resistance_phloem_to_xylem_sh;
  double resistance_phloem_to_xylem_st;
  double lower_bound_W;
  double tau_emergancy;
  double b0_repo;
  double b1_repo;
  double b2_repo;
  double uk_repo;
  double eki_repo;
  double stem_no;
  double m_R_tot;
  double diameter_start_day;
  double GPP_mean;
  double GPP_initial;
};

struct CASSIA_common {
  double a;
  double b;
  double TR0;
  double abs_zero;
  double b_s;
  double theetta_FC;
  double phi_e;
  double K_sat;
  double R_length;
  double M_H20;
  double r_cyl;
  double r_root;
  double ypsilon;
  double Rg_N;
  double Rg_S;
  double Rg_R;
  double gas_const;
  double M_C;
  double M_H;
  double M_O;
  double osmotic_sugar_conc;
  double m_N;
  double Uggla;
};

struct CASSIA_ratios {
  double form_factor;
  double needle_fineroot_ratio;
  double sapwood_share;
  double height_growth_coefficient;
  double diameter_growth_coefficient;
  double height_growth_coefficient_max;
  double height_growth_coefficient_min;
  double diameter_growth_coefficient_max;
  double diameter_growth_coefficient_min;
};

struct phydro_canopy_parameters {
  // Photosynthesis
  double alpha; // Cost for photosynthesis
  double gamma; // Cost for water
  double infra_translation; // Conversion from area biomass ratio to nitrogen price
  double kphio; // Quantum yield
  double rd; // Dark respiration
  double a_jmax; // Nitorgen to jmax ratio

  double p50_leaf;        ///< Leaf or whole-plant hydraulic vulnerability [MPa] (calculated from Xylem P50 and Safety margin)
  double K_leaf;          ///< Leaf conductivity [m]
  double b_leaf;          ///< Shape parameter of leaf vulnerabilty curve [-]
  double cbio;            // kg biomass per mol CO2 = 12.011 gC / mol CO2 * 1e-3 kgC/gC * 2.04 kg biomass/kg

  // Environment
  int n_layers;
  double total_crown_area;
  std::vector<double> z_star;
  std::vector<double> fapar_tot;
  std::vector<double> canopy_openness;

  // Canopy
  double m;
  double n;
  double zm_H; // qm, zm_H Precomputed Geometric parameters
  double qm;
  double fg; // fg upper canopy gap fraction
  double k_light;

  // Weather
  double tau_weather;
  double dt;
};

/*
 * PRELES structure
 */

//struct p1 {
//  double soildepth = 413.0;
//  double ThetaFC = 0.45;
//  double ThetaPWP = 0.118;
//  double tauDrainage = 3.0;
//};

//  GPP model
//struct p2 {
//  double beta = 0.7457;
//  double tau = 10.93;
//  double S0 = -3.063;
//  double Smax = 17.72;
//  double kappa = -0.1027;
//  double gamma = 0.03673;
//  double soilthres = 0.7779; // used for fW with ETmodel = 2 | 4 | 6
//  double bCO2 = 0.5; // used for fW with ETmodel = 1 | 3 | 5
//  double xCO2 = -0.364; // used for fW with ETmodel = 1 | 3 | 5;
//  double t0 = -999; // Birch phenology parameters: 26th Feb = 57 DOY
//  double tcrit = -999; // (Linkosalo et al)           1.5 C
//  double tsumcrit = -999; //                                 134  C
//};

// ET-model
//struct p3 {
//  double beta = 0.2715;
//  double kappa = 0.8351;
//  double chi = 0.07348;
//  double soilthres = 0.9996; // used for fW with ETmodel = 2 | 4
//  double nu = 0.4428;
//};

// Rain and Snow models: interception and melting of snow
//struct p4 {
//  double MeltCoef = 1.2;
  // double Ifrac;
  //  double I0 = 0.33;
  //  double CWmax = 4.970496;
  //  double SnowThreshold = 0.0;
  //  double T_0 = 0.0;
  //};

// Storage components
//struct p5 {
//  double SW = 160.0; // Soilw water at beginning
//  double CW = 0.0; // Canopy water
//  double SOG = 0.0; // Snow on Ground
//  double S = 20.0; // State of temperature acclimation
//};

//struct p6 {
//  double cvGPP; // Coefficients of variation for GPP, ET and SW
//  double cvET;  // Used in MCMC-calibration only
//  double cvSW;
//};

struct p7 {
  double a = 0.1; // gradient
  double b = 0; // intercept
};

#endif
