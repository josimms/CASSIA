#ifndef FINAL_PARAMETERS_H
#define FINAL_PARAMETERS_H

/*
 * CASSIA Structures
 */

struct CASSIA_parameters {
  double Q10_N = 1.898;
  double Rm_N = 0.00267;
  double Q10_S = 1.74788;
  double Rm_S = 0.000055576;
  double Q10_R = 2.5575;
  double Rm_R = 0.00958;
  double sR0 = -2.5;
  double sRc = 30.0;
  double growth_myco = 0.1;
  double root_lifetime = 1.7;
  double HH0 = 10.0;
  double sH0 = -1.359200388;
  double LH = 8.226401284;
  double LH0 = 8.226401284;
  double sHc = 14.596362790;
  double sN0 = -8.37584;
  double LN = 1.849493;
  double LN0 = 1.849493;
  double sNc = 5.263883;
  double HN0 = 1.0;
  double sD0_Trad = -3.724083738;
  double LD;
  double LD0 = 1.293443902;
  double sDc = 5.077004992;
  double sDc_T_count;
  double tau_Ee = 10.686858770;
  double tau_El = 8.789131263;
  double tau_We = 25.294488570;
  double tau_Wl = 35.121486870;
  double tau_GPP = 5.0;
  double Uggla = 1.950000000;
  double sB0 = 171.0;
  double sBc = 85.0;
  double LB = 0.005;
  double cell_d_ew = 0.000035700;
  double cell_d_lw = 0.000024200;
  double cell_l_ew = 0.00259;
  double cell_l_lw = 0.00273;
  double cell_wall_density_ew = 57.0;
  double cell_wall_density_lw = 57.0;
  double wall_thickness_ew = 0.000002610;
  double wall_thickness_lw = 0.000005230;
  double cell_volume_growth_per_day_ew; // TODO: values
  double cell_volume_growth_per_day_lw;
  double density_tree = 400.0;
  double carbon_share = 0.5;
  double D0 = 0.175;
  double h0 = 17.90;
  double n_age = 3.0;
  double n_length = 34.241;
  double h_increment = 309.0938;
  double SLA = 13.0;
  double LR0;
  double mycorrhiza_threshold = 0.025;
  double starch0 = 0.3246781;
  double sugar0 = 0.4184208;
  double starch_needles0 = 0.03;
  double starch_phloem0 = 0.037;
  double starch_xylem_sh0 = 0.034;
  double starch_xylem_st0 = 0.166;
  double starch_roots0 = 0.057;
  double sugar_needles0 = 0.087;
  double sugar_phloem0 = 0.27;
  double sugar_roots0 = 0.014;
  double sugar_xylem_sh0 = 0.0249;
  double sugar_xylem_st0 = 0.021;
  double Wala_needles = 0.0;
  double Wala_phloem = 0.0;
  double Wala_xylem_sh = 0.0;
  double Wala_xylem_st = 0.0;
  double Wala_roots = 0.0;
  double carbon_sugar = 0.4211;
  double carbon_starch = 0.4444;
  double alfa_needles = 3.0;
  double alfa_phloem = 3.0;
  double alfa_xylem_sh = 3.0;
  double alfa_xylem_st = 3.0;
  double alfa_roots = 3.0;
  double tau_s = 2.0;
  double tau_t = 2.0;
  double starch00 = 0.3246781;
  double sugar00 = 0.4184208;
  double Wala = 0.0;
  double alfa = 3.0;
  double Q10s = 3.0;
  double Q10d = 1.8;
  double SCb = 0.23;
  double sugar_level = 0.41;
  double Ad0_needles = 0.017;
  double Ad0_phloem = 0.008;
  double Ad0_roots = 0.0002;
  double Ad0_xylem_sh = 0.0002;
  double Ad0_xylem_st = 0.047;
  double lambda_needles = 0.197;
  double lambda_phloem = 0.05301;
  double lambda_roots = 0.211;
  double lambda_xylem_sh = 0.00401;
  double lambda_xylem_st = 0.00401;
  double delta_needles = 0.729;
  double delta_phloem = 0.832;
  double delta_roots = 0.853;
  double delta_xylem_sh = 0.762;
  double delta_xylem_st = 0.294;
  double lower_bound_needles = 0.02;
  double lower_bound_phloem = 0.03;
  double lower_bound_roots = 0.05;
  double lower_bound_xylem_sh = 0.03;
  double lower_bound_xylem_st = 0.1;
  double tau_emergancy_needles = 3;
  double tau_emergancy_phloem = 3;
  double tau_emergancy_roots = 3;
  double tau_emergancy_xylem_sh = 3;
  double tau_emergancy_xylem_st = 3;
  double resistance_needles_to_phloem = 0.3;
  double resistance_phloem_to_roots = 0.072;
  double resistance_phloem_to_xylem_sh = 0.188;
  double resistance_phloem_to_xylem_st = 0.17;
  double lower_bound_W;
  double tau_emergancy;
  double b0_repo = -6.303;
  double b1_repo = -14.472;
  double b2_repo = -3.976;
  double uk_repo = 0.109;
  double eki_repo = -0.118;
  double stem_no = 1010;
  double m_R_tot = 0; // TODO: This shouldn't really exist here, not quite sure here
};

struct CASSIA_common {
  double a = 0.185;
  double b = 18.4;
  double TR0 = 0;
  double abs_zero = 273.15;
  double b_s = 4.14;
  double theetta_FC = 0.62;
  double phi_e = -0.00068;
  double K_sat = 24.5;
  double R_length = 5300;
  double M_H20 = 0.018;
  double r_cyl = 0.00425;
  double r_root = 0.003;
  double ypsilon = 1e-14;
  double Rg_N = 0.35;
  double Rg_S = 0.3;
  double Rg_R = 0.35;
  double gas_const = 0.314;
  double M_C = 12.01;
  double M_H = 1.008;
  double M_O = 16;
  double osmotic_sugar_conc = 2e+06;
  double m_N;
  double Uggla;
};

struct CASSIA_ratios {
  double form_factor = 0.6;
  double needle_fineroot_ratio;
  double sapwood_share = 0.8;
  double height_growth_coefficient = 4.3;
  double diameter_growth_coefficient = 1.6;
  double height_growth_coefficient_max = 5.5;
  double height_growth_coefficient_min = 3.8;
  double diameter_growth_coefficient_max = 1.9;
  double diameter_growth_coefficient_min = 1.5;
};

/*
 * PRELES structure
 */

struct p1 {
  double soildepth = 413.0;
  double ThetaFC = 0.45;
  double ThetaPWP = 0.118;
  double tauDrainage = 3.0;
};

//  GPP model
struct p2 {
  double beta = 0.7457;
  double tau = 10.93;
  double S0 = -3.063;
  double Smax = 17.72;
  double kappa = -0.1027;
  double gamma = 0.03673;
  double soilthres = 0.7779; // used for fW with ETmodel = 2 | 4 | 6
  double bCO2 = 0.5; // used for fW with ETmodel = 1 | 3 | 5
  double xCO2 = -0.364; // used for fW with ETmodel = 1 | 3 | 5;
  double t0 = -999; // Birch phenology parameters: 26th Feb = 57 DOY
  double tcrit = -999; // (Linkosalo et al)           1.5 C
  double tsumcrit = -999; //                                 134  C
};

// ET-model
struct p3 {
  double beta = 0.2715;
  double kappa = 0.8351;
  double chi = 0.07348;
  double soilthres = 0.9996; // used for fW with ETmodel = 2 | 4
  double nu = 0.4428;
};

// Rain and Snow models: interception and melting of snow
struct p4 {
  double MeltCoef = 1.2;
  // double Ifrac;
  double I0 = 0.33;
  double CWmax = 4.970496;
  double SnowThreshold = 0.0;
  double T_0 = 0.0;
};

// Storage components
struct p5 {
  double SW = 160.0; // Soilw water at beginning
  double CW = 0.0; // Canopy water
  double SOG = 0.0; // Snow on Ground
  double S = 20.0; // State of temperature acclimation
};

struct p6 {
  double cvGPP; // Coefficients of variation for GPP, ET and SW
  double cvET;  // Used in MCMC-calibration only
  double cvSW;
};

struct p7 {
  double a = 0.1; // gradient
  double b = 0; // intercept
};

#endif
