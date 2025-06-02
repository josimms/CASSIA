#include "CASSIA.h"

// [[Rcpp::export]]
void printColumnNames(Rcpp::DataFrame df) {
  if (df.size() == 0) {
    Rcpp::Rcout << "The DataFrame is empty." << std::endl;
    return;
  }

  Rcpp::CharacterVector col_names = df.names();

  Rcpp::Rcout << "Column names:" << std::endl;
  for (int i = 0; i < col_names.length(); ++i) {
    Rcpp::Rcout << col_names[i] << std::endl;
  }
}

CASSIA_parameters make_CASSIA_parameters(Rcpp::DataFrame input_parameters,
                                         Rcpp::DataFrame input_sperling) {
  CASSIA_parameters out;
  std::vector<double> temp = input_parameters["Q10.N"];
  out.Q10_N = temp[0];
  std::vector<double> temp1 = input_parameters["Rm.NR"];
  out.Rm_N = temp1[0];
  std::vector<double> temp2 = input_parameters["Q10.S"];
  out.Q10_S = temp2[0];
  std::vector<double> temp3 = input_parameters["Rm.S"];
  out.Rm_S = temp3[0];
  std::vector<double> temp4 = input_parameters["Q10.R"];
  out.Q10_R = temp4[0];
  std::vector<double> temp5 = input_parameters["Rm.R"];
  out.Rm_R = temp5[0];
  std::vector<double> temp6 = input_parameters["sR0"];
  out.sR0 = temp6[0];
  std::vector<double> temp7 = input_parameters["sRc"];
  out.sRc = temp7[0];
  std::vector<double> temp8 = input_parameters["growth.myco"];
  out.growth_myco = temp8[0];
  std::vector<double> temp9 = input_parameters["root.lifetime"];
  out.root_lifetime = temp9[0];
  std::vector<double> temp10 = input_parameters["HH0"];
  out.HH0 = temp10[0];
  std::vector<double> temp11 = input_parameters["sH0"];
  out.sH0 = temp11[0];
  std::vector<double> temp12 = input_parameters["LH"];
  out.LH = temp12[0];
  std::vector<double> temp13 = input_parameters["LH0"];
  out.LH0 = temp13[0];
  std::vector<double> temp14 = input_parameters["sHc"];
  out.sHc = temp14[0];
  std::vector<double> temp15 = input_parameters["sN0"];
  out.sN0 = temp15[0];
  std::vector<double> temp16 = input_parameters["LN"];
  out.LN = temp16[0];
  std::vector<double> temp17 = input_parameters["LN0"];
  out.LN0 = temp17[0];
  std::vector<double> temp18 = input_parameters["sNc"];
  out.sNc = temp18[0];
  std::vector<double> temp19 = input_parameters["HN0"];
  out.HN0 = temp19[0];
  std::vector<double> temp20 = input_parameters["sD0.Trad"];
  out.sD0_Trad = temp20[0];
  std::vector<double> temp21 = input_parameters["LD"];
  out.LD = temp21[0];
  std::vector<double> temp22 = input_parameters["LD0"];
  out.LD0 = temp22[0];
  std::vector<double> temp23 = input_parameters["sDc"];
  out.sDc = temp23[0];
  std::vector<double> temp24 = input_parameters["sDc.T.count"];
  out.sDc_T_count = temp24[0];
  std::vector<double> temp25 = input_parameters["tau.Ee"];
  out.tau_Ee = temp25[0];
  std::vector<double> temp26 = input_parameters["tau.El"];
  out.tau_El = temp26[0];
  std::vector<double> temp27 = input_parameters["tau.We"];
  out.tau_We = temp27[0];
  std::vector<double> temp28 = input_parameters["tau.Wl"];
  out.tau_Wl = temp28[0];
  std::vector<double> temp29 = input_parameters["tau.GPP"];
  out.tau_GPP = temp29[0];
  std::vector<double> temp30 = input_parameters["Uggla"];
  out.Uggla = temp30[0];
  std::vector<double> temp31 = input_parameters["sB0"];
  out.sB0 = temp31[0];
  std::vector<double> temp32 = input_parameters["sBc"];
  out.sBc = temp32[0];
  std::vector<double> temp33 = input_parameters["LB"];
  out.LB = temp33[0];
  std::vector<double> temp34 = input_parameters["cell.d.ew"];
  out.cell_d_ew = temp34[0];
  std::vector<double> temp35 = input_parameters["cell.d.lw"];
  out.cell_d_lw = temp35[0];
  std::vector<double> temp36 = input_parameters["cell.l.ew"];
  out.cell_l_ew = temp36[0];
  std::vector<double> temp37 = input_parameters["cell.l.lw"];
  out.cell_l_lw = temp37[0];
  std::vector<double> temp38 = input_parameters["cell.wall.density.ew"];
  out.cell_wall_density_ew = temp38[0];
  std::vector<double> temp39 = input_parameters["cell.wall.density.lw"];
  out.cell_wall_density_lw = temp39[0];
  std::vector<double> temp40 = input_parameters["wall.thickness.ew"];
  out.wall_thickness_ew = temp40[0];
  std::vector<double> temp41 = input_parameters["wall.thickness.lw"];
  out.wall_thickness_lw = temp41[0];
  std::vector<double> temp42 = input_parameters["cell.volume.growth.per.day.ew"]; // TODO: values
  out.cell_volume_growth_per_day_ew = temp42[0];
  std::vector<double> temp43 = input_parameters["cell.volumn.growth.per.day.lw"];
  out.cell_volume_growth_per_day_lw = temp43[0];
  std::vector<double> temp44 = input_parameters["density_tree"];
  out.density_tree = temp44[0];
  std::vector<double> temp45 = input_parameters["carbon_share"];
  out.carbon_share = temp45[0];
  std::vector<double> temp46 = input_parameters["D0"];
  out.D0 = temp46[0];
  std::vector<double> temp47 = input_parameters["h0"];
  out.h0 = temp47[0];
  std::vector<double> temp48 = input_parameters["n_age"];
  out.n_age = temp48[0];
  std::vector<double> temp49 = input_parameters["n_lenght"];
  out.n_length = temp49[0];
  std::vector<double> temp50 = input_parameters["h_increment"];
  out.h_increment = temp50[0];
  std::vector<double> temp51 = input_parameters["SLA"];
  out.SLA = temp51[0];
  std::vector<double> temp52 = input_parameters["LR0"];
  out.LR0 = temp52[0];
  std::vector<double> temp53 = input_sperling["myco.thresh"];
  out.mycorrhiza_threshold = temp53[0];
  std::vector<double> temp54 = input_sperling["starch0"];
  out.starch0 = temp54[0];
  std::vector<double> temp55 = input_sperling["sugar0"];
  out.sugar0 = temp55[0];
  std::vector<double> temp56 = input_sperling["starch.needles0"];
  out.starch_needles0 = temp56[0];
  std::vector<double> temp57 = input_sperling["starch.phloem0"];
  out.starch_phloem0 = temp57[0];
  std::vector<double> temp58 = input_sperling["starch.xylem.sh0"];
  out.starch_xylem_sh0 = temp58[0];
  std::vector<double> temp59 = input_sperling["starch.xylem.st0"];
  out.starch_xylem_st0 = temp59[0];
  std::vector<double> temp60 = input_sperling["starch.roots0"];
  out.starch_roots0 = temp60[0];
  std::vector<double> temp61 = input_sperling["sugar.needles0"];
  out.sugar_needles0 = temp61[0];
  std::vector<double> temp62 = input_sperling["sugar.phloem0"];
  out.sugar_phloem0 = temp62[0];
  std::vector<double> temp63 = input_sperling["sugar.roots0"];
  out.sugar_roots0 = temp63[0];
  std::vector<double> temp64 = input_sperling["sugar.xylem.sh0"];
  out.sugar_xylem_sh0 = temp64[0];
  std::vector<double> temp65 = input_sperling["sugar.xylem.st0"];
  out.sugar_xylem_st0 = temp65[0];
  std::vector<double> temp66 = input_sperling["Wala.needles"];
  out.Wala_needles = temp66[0];
  std::vector<double> temp67 = input_sperling["Wala.phloem"];
  out.Wala_phloem = temp67[0];
  std::vector<double> temp68 = input_sperling["Wala.xylem.sh"];
  out.Wala_xylem_sh = temp68[0];
  std::vector<double> temp69 = input_sperling["Wala.xylem.st"];
  out.Wala_xylem_st = temp69[0];
  std::vector<double> temp70 = input_sperling["Wala.roots"];
  out.Wala_roots = temp70[0];
  std::vector<double> temp71 = input_sperling["carbon.sugar"];
  out.carbon_sugar = temp71[0];
  std::vector<double> temp72 = input_sperling["carbon.starch"];
  out.carbon_starch = temp72[0];
  std::vector<double> temp73 = input_sperling["alfa.needles"];
  out.alfa_needles = temp73[0];
  std::vector<double> temp74 = input_sperling["alfa.phloem"];
  out.alfa_phloem = temp74[0];
  std::vector<double> temp75 = input_sperling["alfa.xylem.sh"];
  out.alfa_xylem_sh = temp75[0];
  std::vector<double> temp76 = input_sperling["alfa.xylem.st"];
  out.alfa_xylem_st = temp76[0];
  std::vector<double> temp77 = input_sperling["alfa.roots"];
  out.alfa_roots = temp77[0];
  std::vector<double> temp78 = input_sperling["tau.s"];
  out.tau_s = temp78[0];
  std::vector<double> temp79 = input_sperling["tau.t"];
  out.tau_t = temp79[0];
  std::vector<double> temp80 = input_sperling["starch00"];
  out.starch00 = temp80[0];
  std::vector<double> temp81 = input_sperling["sugar00"];
  out.sugar00 = temp81[0];
  std::vector<double> temp82 = input_sperling["Wala"];
  out.Wala = temp82[0];
  std::vector<double> temp83 = input_sperling["alfa"];
  out.alfa = temp83[0];
  std::vector<double> temp84 = input_sperling["Q10s"];
  out.Q10s = temp84[0];
  std::vector<double> temp85 = input_sperling["Q10d"];
  out.Q10d = temp85[0];
  std::vector<double> temp86 = input_sperling["SCb"];
  out.SCb = temp86[0];
  std::vector<double> temp87 = input_sperling["sugar.level"];
  out.sugar_level = temp87[0];
  std::vector<double> temp88 = input_sperling["Ad0.needles"];
  out.Ad0_needles = temp88[0];
  std::vector<double> temp89 = input_sperling["Ad0.phloem"];
  out.Ad0_phloem = temp89[0];
  std::vector<double> temp90 = input_sperling["Ad0.roots"];
  out.Ad0_roots = temp90[0];
  std::vector<double> temp91 = input_sperling["Ad0.xylem.sh"];
  out.Ad0_xylem_sh = temp91[0];
  std::vector<double> temp92 = input_sperling["Ad0.xylem.st"];
  out.Ad0_xylem_st = temp92[0];
  std::vector<double> temp93 = input_sperling["lamda.needles"];
  out.lambda_needles = temp93[0];
  std::vector<double> temp94 = input_sperling["lamda.phloem"];
  out.lambda_phloem = temp94[0];
  std::vector<double> temp95 = input_sperling["lamda.roots"];
  out.lambda_roots = temp95[0];
  std::vector<double> temp96 = input_sperling["lamda.xylem.sh"];
  out.lambda_xylem_sh = temp96[0];
  std::vector<double> temp97 = input_sperling["lamda.xylem.st"];
  out.lambda_xylem_st = temp97[0];
  std::vector<double> temp98 = input_sperling["delta.needles"];
  out.delta_needles = temp98[0];
  std::vector<double> temp99 = input_sperling["delta.phloem"];
  out.delta_phloem = temp99[0];
  std::vector<double> temp100 = input_sperling["delta.roots"];
  out.delta_roots = temp100[0];
  std::vector<double> temp101 = input_sperling["delta.xylem.sh"];
  out.delta_xylem_sh = temp101[0];
  std::vector<double> temp102 = input_sperling["delta.xylem.st"];
  out.delta_xylem_st = temp102[0];
  std::vector<double> temp110 = input_parameters["lower_bound_needles"];
  out.lower_bound_needles = temp110[0];
  std::vector<double> temp111 = input_parameters["lower_bound_phloem"];
  out.lower_bound_phloem = temp111[0];
  std::vector<double> temp112 = input_parameters["lower_bound_roots"];
  out.lower_bound_roots = temp112[0];
  std::vector<double> temp113 = input_parameters["lower_bound_xylem_sh"];
  out.lower_bound_xylem_sh = temp113[0];
  std::vector<double> temp114 = input_parameters["lower_bound_xylem_st"];
  out.lower_bound_xylem_st = temp114[0];
  std::vector<double> temp115 = input_parameters["tau_emergancy_needles"];
  out.tau_emergancy_needles = temp115[0];
  std::vector<double> temp116 = input_parameters["tau_emergancy_phloem"];
  out.tau_emergancy_phloem = temp116[0];
  std::vector<double> temp117 = input_parameters["tau_emergancy_roots"];
  out.tau_emergancy_roots = temp117[0];
  std::vector<double> temp118 = input_parameters["tau_emergancy_xylem_sh"];
  out.tau_emergancy_xylem_sh = temp118[0];
  std::vector<double> temp119 = input_parameters["tau_emergancy_xylem_st"];
  out.tau_emergancy_xylem_st = temp119[0];
  std::vector<double> temp103 = input_sperling["percentage_needle_storage"];
  out.percentage_needle_storage = temp103[0];
  std::vector<double> temp104 = input_sperling["percentage_xylem_sh_storage"];
  out.percentage_xylem_sh_storage = temp104[0];
  std::vector<double> temp105 = input_sperling["percentage_xylem_st_storage"];
  out.percentage_xylem_st_storage = temp105[0];
  std::vector<double> temp106 = input_sperling["percentage_phloem_storage"];
  out.percentage_phloem_storage = temp106[0];
  std::vector<double> temp1065 = input_sperling["percentage_roots_storage"];
  out.percentage_roots_storage = temp1065[0];
  std::vector<double> temp120 = input_parameters["lower_bound_W"];
  out.lower_bound_W = temp120[0];
  std::vector<double> temp121 = input_parameters["tau_emergancy"];
  out.tau_emergancy = temp121[0];
  std::vector<double> temp107 = input_parameters["b0_repo"];
  out.b0_repo = temp107[0];
  std::vector<double> temp108 = input_parameters["b1_repo"];
  out.b1_repo = temp108[0];
  std::vector<double> temp109 = input_parameters["b2_repo"];
  out.b2_repo = temp109[0];
  std::vector<double> temp122 = input_parameters["uk_repo"];
  out.uk_repo = temp122[0];
  std::vector<double> temp123 = input_parameters["eki_repo"];
  out.eki_repo = temp123[0];
  std::vector<double> temp124 = input_parameters["stem_no"];
  out.stem_no = temp124[0];
  std::vector<double> temp125 = input_parameters["diameter_start_day"];
  out.diameter_start_day = temp125[0];
  std::vector<double> temp126 = input_parameters["GPP_mean"];
  out.GPP_mean = temp126[0];
  std::vector<double> temp127 = input_parameters["GPP_initial"];
  out.GPP_initial = temp127[0];
  return out;
}

CASSIA_common make_common(Rcpp::DataFrame input) {
  CASSIA_common out;
  std::vector<double> temp = input["a"];
  out.a = temp[0];
  std::vector<double> temp1 = input["b"];
  out.b = temp1[0];
  std::vector<double> temp2 = input["TR0"];
  out.TR0 = temp2[0];
  std::vector<double> temp3 = input["abs_zero"];
  out.abs_zero = temp3[0];
  std::vector<double> temp4 = input["b.s"];
  out.b_s = temp4[0];
  std::vector<double> temp5 = input["theetta.FC"];
  out.theetta_FC = temp5[0];
  std::vector<double> temp6 = input["phi.e"];
  out.phi_e = temp6[0];
  std::vector<double> temp7 = input["K.sat"];
  out.K_sat = temp7[0];
  std::vector<double> temp8 = input["R.length"];
  out.R_length = temp8[0];
  std::vector<double> temp9 = input["M.H20"];
  out.M_H20 = temp9[0];
  std::vector<double> temp10 = input["r.cyl"];
  out.r_cyl = temp10[0];
  std::vector<double> temp11 = input["r.root"];
  out.r_root = temp11[0];
  std::vector<double> temp12 = input["ypsilon"];
  out.ypsilon = temp12[0];
  std::vector<double> temp13 = input["Rg.N"];
  out.Rg_N = temp13[0];
  std::vector<double> temp14 = input["Rg.S"];
  out.Rg_S = temp14[0];
  std::vector<double> temp15 = input["Rg.R"];
  out.Rg_R = temp15[0];
  std::vector<double> temp16 = input["gas.const"];
  out.gas_const = temp16[0];
  std::vector<double> temp17 = input["M.C"];
  out.M_C = temp17[0];
  std::vector<double> temp18 = input["M.H"];
  out.M_H = temp18[0];
  std::vector<double> temp19 = input["M.O"];
  out.M_O = temp19[0];
  std::vector<double> temp20 = input["osmotic.sugar.conc"];
  out.osmotic_sugar_conc = temp20[0];
  std::vector<double> temp22 = input["m_N"];
  out.m_N = temp22[0];
  std::vector<double> temp23 = input["Uggla"];
  out.Uggla = temp23[0];
  return(out);
}


CASSIA_common read_common_parameters(const std::string& filename) {
  CASSIA_common common_params;

  // Map to store parameter names and values
  std::unordered_map<std::string, double> params;

  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return common_params;
  }

  std::string line;
  std::getline(file, line); // Read header line

  // Read each line of the file
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string parameter;
    std::string value_str;

    // Read parameter and value from the line
    std::getline(ss, parameter, ',');
    std::getline(ss, value_str);

    // Convert value to double and store in the map
    double value = std::stod(value_str);
    params[parameter] = value;
  }
  file.close();

  // Assign values to the structure based on the map
  common_params.a = params["a"];
  common_params.b = params["b"];
  common_params.TR0 = params["TR0"];
  common_params.abs_zero = params["abs_zero"];
  common_params.b_s = params["b_s"];
  common_params.theetta_FC = params["theetta_FC"];
  common_params.phi_e = params["phi_e"];
  common_params.K_sat = params["K_sat"];
  common_params.R_length = params["R_length"];
  common_params.M_H20 = params["M_H20"];
  common_params.r_cyl = params["r_cyl"];
  common_params.r_root = params["r_root"];
  common_params.ypsilon = params["ypsilon"];
  common_params.Rg_N = params["Rg_N"];
  common_params.Rg_S = params["Rg_S"];
  common_params.Rg_R = params["Rg_R"];
  common_params.gas_const = params["gas_const"];
  common_params.M_C = params["M_C"];
  common_params.M_H = params["M_H"];
  common_params.M_O = params["M_O"];
  common_params.osmotic_sugar_conc = params["osmotic_sugar_conc"];
  common_params.m_N = params["m_N"];
  common_params.Uggla = params["Uggla"];

  return common_params;
}

CASSIA_ratios read_ratios(const std::string& filename, const std::string& site) {
  CASSIA_ratios ratios;
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return ratios;
  }

  std::string line, header;
  std::getline(file, header);  // Skip the header line

  // Determine the column index for the chosen site
  std::map<std::string, int> site_index = {
    {"Hyde", 1},
    {"Lettosuo", 2},
    {"china", 3},
    {"vario", 4}
  };

  int column_index = site_index[site];

  while (std::getline(file, line)) {
    std::istringstream linestream(line);
    std::string parameter, value_str;
    double value;

    // Read the parameter name
    std::getline(linestream, parameter, ',');

    // Skip columns until we reach the desired site
    for (int i = 0; i < column_index; ++i) {
      std::getline(linestream, value_str, ',');
    }

    // Convert the value string to a double
    if (!value_str.empty()) {
      value = std::stod(value_str);
    } else {
      value = 0.0; // Handle missing values if necessary
    }

    // Assign values to the corresponding fields in the structure
    if (parameter == "form_factor") {
      ratios.form_factor = value;
    } else if (parameter == "needle_fineroot_ratio") {
      ratios.needle_fineroot_ratio = value;
    } else if (parameter == "sapwood.share" || parameter == "sapwood_share") {
      ratios.sapwood_share = value;
    } else if (parameter == "height_growth_coefficient") {
      ratios.height_growth_coefficient = value;
    } else if (parameter == "diameter_growth_coefficient") {
      ratios.diameter_growth_coefficient = value;
    } else if (parameter == "height_growth_coefficient_max") {
      ratios.height_growth_coefficient_max = value;
    } else if (parameter == "height_growth_coefficient_min") {
      ratios.height_growth_coefficient_min = value;
    } else if (parameter == "diameter_growth_coefficient_max") {
      ratios.diameter_growth_coefficient_max = value;
    } else if (parameter == "diameter_growth_coefficient_min") {
      ratios.diameter_growth_coefficient_min = value;
    }
  }

  file.close();
  return ratios;
}


CASSIA_ratios make_ratios(Rcpp::DataFrame input) {
  CASSIA_ratios out;
  std::vector<double> temp = input["form_factor"];
  out.form_factor = temp[0];
  std::vector<double> temp1 = input["needle_fineroot_ratio"];
  out.needle_fineroot_ratio = temp1[0];
  std::vector<double> temp2 = input["sapwood.share"];
  out.sapwood_share = temp2[0];
  std::vector<double> temp3 = input["height_growth_coefficient"];
  out.height_growth_coefficient = temp3[0];
  std::vector<double> temp4 = input["diameter_growth_coefficient"];
  out.diameter_growth_coefficient = temp4[0];
  std::vector<double> temp5 = input["height_growth_coefficient_max"];
  out.height_growth_coefficient_max = temp5[0];
  std::vector<double> temp6 = input["height_growth_coefficient_min"];
  out.height_growth_coefficient_min = temp6[0];
  std::vector<double> temp7 = input["diameter_growth_coefficient_max"];
  out.diameter_growth_coefficient_max = temp7[0];
  std::vector<double> temp8 = input["diameter_growth_coefficient_min"];
  out.diameter_growth_coefficient_min = temp8[0];
  return(out);
}

p1 make_p1(std::vector<double> input) {
  p1 out;
  out.soildepth = input[0];
  out.ThetaFC = input[1];
  out.ThetaPWP = input[2];
  out.tauDrainage = input[3];
  return(out);
}

p2 make_p2(std::vector<double> input) {
  p2 out;
  out.beta = input[4];
  out.tau = input[5];
  out.S0 = input[6];
  out.Smax = input[7];
  out.kappa = input[8];
  out.gamma = input[9];
  out.soilthres = input[10];
  out.bCO2 = input[11];
  out.xCO2 = input[12];
  out.t0 = input[27];
  out.tcrit = input[28];
  out.tsumcrit = input[29];
  return(out);
}

p3 make_p3(std::vector<double> input) {
  p3 out;
  out.beta = input[13];
  out.kappa = input[14];
  out.chi = input[15];
  out.soilthres = input[16];
  out.nu = input[17];
  return(out);
}

p4 make_p4(std::vector<double> input) {
  p4 out;
  out.MeltCoef = input[18];
  out.I0 = input[19];
  out.CWmax = input[20];
  out.SnowThreshold = input[21];
  out.T_0 = input[22];
  return(out);
}

p5 make_p5(std::vector<double> input) {
  p5 out;
  out.SW = input[23];
  out.CW = input[24];
  out.SOG = input[25];
  out.S = input[26];
  return(out);
}

p7 make_p7(std::vector<double> input) {
  p7 out;
  out.a = input[30];
  out.b = input[31];
  return(out);
}

phydro_canopy_parameters parPhydro_initalise(std::vector<double> phydro_params) {
  phydro_canopy_parameters parPhydro;
  parPhydro.alpha = phydro_params[0];           // Cost of maintaining photosynthetic capacity (Ref: Joshi et al 2022, removed shrubs and gymnosperms)
  parPhydro.gamma = phydro_params[1];           // Cost of maintaining hydraulic pathway  (Ref: Joshi Slack)
  parPhydro.infra_translation = phydro_params[2];    // Translation from biomass area ratio to Ib
  parPhydro.kphio = phydro_params[3]; // 0.087            // Quantum yield efficiency
  parPhydro.rd = phydro_params[4]; // 0.015              // ratio of leaf dark respiration rate to vcmax [-]  (Ref: 0.011 in Farquhar et al 1980, 0.015 in Collatz et al 1991)
  parPhydro.a_jmax = phydro_params[5]; // TODO: 800 or 50?

  parPhydro.k_light = phydro_params[6];

  // Note the p50_leaf was commented and this was the value p50_xylem -2.29 was there
  parPhydro.p50_leaf = phydro_params[7];          // Leaf P50 [MPa]
  parPhydro.K_leaf = phydro_params[8];         // Leaf conductance [m]  ---> ** Calibrated to gs **
  parPhydro.b_leaf = phydro_params[9];               // Shape parameter of xylem vulnerabilty curve [-]
  parPhydro.cbio = phydro_params[10];           // kg biomass per mol CO2 = 12.011 gC / mol CO2 * 1e-3 kgC/gC * 2.04 kg biomass/kg
  // TODO: replace with a boreal shape
  parPhydro.m = phydro_params[11];                  // crown shape smoothness
  parPhydro.n = phydro_params[12];                    // crown top-heaviness
  parPhydro.fg = phydro_params[13];                 // upper canopy gap fraction
  // qm defined in plant_architecture
  parPhydro.qm = parPhydro.m * parPhydro.n * pow((parPhydro.n - 1) / (parPhydro.m * parPhydro.n - 1), 1 - 1 / parPhydro.n) * pow((parPhydro.m - 1) * parPhydro.n / (parPhydro.m * parPhydro.n - 1), parPhydro.m - 1);
  // zm_H defined in plant_architecture
  parPhydro.zm_H = pow((parPhydro.n - 1) / (parPhydro.m * parPhydro.n - 1), 1 / parPhydro.n);

  // This value is random, although 15 is one of the values from PlantFate
  parPhydro.z_star = {phydro_params[14], phydro_params[15], phydro_params[16], phydro_params[17]};
  parPhydro.canopy_openness = {phydro_params[18], phydro_params[19], phydro_params[20], phydro_params[21]};

  return(parPhydro);
}
// Helper function to convert double to string with specified precision
std::string doubleToStringPrecision(double value, int precision = 6) {
  std::ostringstream out;
  out.precision(precision);
  out << std::fixed << value;
  return out.str();
}

void print_phydro_parameters(const phydro_canopy_parameters& params) {
  // Photosynthesis parameters
  std::cout << "Photosynthesis Parameters:\n";
  std::cout << "  alpha: " << doubleToStringPrecision(params.alpha) << " (Cost for photosynthesis)\n";
  std::cout << "  gamma: " << doubleToStringPrecision(params.gamma) << " (Cost for water)\n";
  std::cout << "  infra_translation: " << doubleToStringPrecision(params.infra_translation) << " (Conversion from area biomass ratio to nitrogen price)\n";
  std::cout << "  kphio: " << doubleToStringPrecision(params.kphio) << " (Quantum yield)\n";
  std::cout << "  rd: " << doubleToStringPrecision(params.rd) << " (Dark respiration)\n";
  std::cout << "  a_jmax: " << doubleToStringPrecision(params.a_jmax) << " (Nitrogen to jmax ratio)\n";
  std::cout << "  p50_leaf: " << doubleToStringPrecision(params.p50_leaf) << " (Leaf hydraulic vulnerability [MPa])\n";
  std::cout << "  K_leaf: " << doubleToStringPrecision(params.K_leaf) << " (Leaf conductivity [m])\n";
  std::cout << "  b_leaf: " << doubleToStringPrecision(params.b_leaf) << " (Shape parameter of leaf vulnerability curve)\n";
  std::cout << "  cbio: " << doubleToStringPrecision(params.cbio) << " (kg biomass per mol CO2)\n";

  // Environment parameters
  std::cout << "\nEnvironment Parameters:\n";
  std::cout << "  n_layers: " << params.n_layers << " (Number of layers)\n";
  std::cout << "  total_crown_area: " << doubleToStringPrecision(params.total_crown_area) << "\n";

  // z_star vector
  std::cout << "  z_star: [";
  for (size_t i = 0; i < params.z_star.size(); ++i) {
    std::cout << doubleToStringPrecision(params.z_star[i]);
    if (i < params.z_star.size() - 1) std::cout << ", ";
  }
  std::cout << "]\n";

  // fapar_tot vector
  std::cout << "  fapar_tot: [";
  for (size_t i = 0; i < params.fapar_tot.size(); ++i) {
    std::cout << doubleToStringPrecision(params.fapar_tot[i]);
    if (i < params.fapar_tot.size() - 1) std::cout << ", ";
  }
  std::cout << "]\n";

  // canopy_openness vector
  std::cout << "  canopy_openness: [";
  for (size_t i = 0; i < params.canopy_openness.size(); ++i) {
    std::cout << doubleToStringPrecision(params.canopy_openness[i]);
    if (i < params.canopy_openness.size() - 1) std::cout << ", ";
  }
  std::cout << "]\n";

  // Canopy parameters
  std::cout << "\nCanopy Parameters:\n";
  std::cout << "  m: " << doubleToStringPrecision(params.m) << " (Canopy shape parameter)\n";
  std::cout << "  n: " << doubleToStringPrecision(params.n) << " (Canopy shape parameter)\n";
  std::cout << "  zm_H: " << doubleToStringPrecision(params.zm_H) << " (Precomputed geometric parameter)\n";
  std::cout << "  qm: " << doubleToStringPrecision(params.qm) << " (Precomputed geometric parameter)\n";
  std::cout << "  fg: " << doubleToStringPrecision(params.fg) << " (Upper canopy gap fraction)\n";
  std::cout << "  k_light: " << doubleToStringPrecision(params.k_light) << " (Light extinction coefficient)\n";

  // Weather parameters
  std::cout << "\nWeather Parameters:\n";
  std::cout << "  tau_weather: " << doubleToStringPrecision(params.tau_weather) << "\n";
  std::cout << "  dt: " << doubleToStringPrecision(params.dt) << " (Time step)\n";
}


parameters_soil parameters_initalise_test(std::vector<double> parameters_R) {
  parameters_soil out;
  out.microbe_turnover = parameters_R[0];
  out.NC_in_root_opt = parameters_R[1];
  out.NC_fungal_opt = parameters_R[2];
  out.NC_microbe_opt = parameters_R[3];
  out.percentage_C_biomass = parameters_R[4];
  // The vectors here are for the three types of nitrogen
  out.N_limits_myco = {parameters_R[5], parameters_R[6], parameters_R[7]};
  out.N_k_myco = {parameters_R[8], parameters_R[9], parameters_R[10]};
  out.SWC_limits_myco = {parameters_R[11], parameters_R[12], parameters_R[13]};
  out.N_limits_plant = {parameters_R[14], parameters_R[15], parameters_R[16]};
  out.N_k_plant = {parameters_R[17], parameters_R[18], parameters_R[19]};
  out.SWC_limits_plant = {parameters_R[20], parameters_R[21], parameters_R[22]};
  out.N_limits_microbes = {parameters_R[23], parameters_R[24], parameters_R[25], parameters_R[26], parameters_R[27]};
  out.N_k_microbes = {parameters_R[28], parameters_R[29], parameters_R[30], parameters_R[31], parameters_R[32]};
  out.SWC_limits_microbes = {parameters_R[33], parameters_R[34], parameters_R[35], parameters_R[36], parameters_R[37]};
  out.C_limits = {parameters_R[38], parameters_R[39], parameters_R[40]};
  out.NH4_on_NO3 = parameters_R[41];
  out.optimal_root_fungal_biomass_ratio = parameters_R[42];
  out.turnover_mantle = parameters_R[43];
  out.turnover_ERM = parameters_R[44];
  out.turnover_roots = parameters_R[45];
  out.turnover_roots_mycorrhized = parameters_R[46];
  out.turnover_fungal = parameters_R[47];
  out.mantle_mass = parameters_R[48];
  out.ERM_mass = parameters_R[49];
  out.growth_C = parameters_R[50];
  out.growth_N = parameters_R[51];
  out.C_value_param_myco = parameters_R[52];
  out.N_value_param_myco = parameters_R[53];
  out.C_value_param_plant = parameters_R[54];
  out.N_value_param_plant = parameters_R[55];
  return(out);
};


MYCOFON_function_out MYCOFON_structure_conversion(Rcpp::List input) {
  MYCOFON_function_out out;
  out.C_biomass = input[0];
  out.C_roots = input[1];
  out.C_fungal = input[2];
  out.N_roots = input[3];
  out.N_fungal = input[4];
  out.uptake_plant = input[5];
  out.uptake_NH4_plant = input[6];
  out.uptake_NO3_plant = input[7];
  out.uptake_Norg_plant = input[8];
  out.uptake_fungal = input[9];
  out.uptake_NH4_fungal = input[10];
  out.uptake_NO3_fungal = input[11];
  out.uptake_Norg_fungal = input[12];
  out.from_CASSIA = input[13];
  out.to_CASSIA = input[14];
  out.Plant_demand = input[15];
  out.Fungal_demand = input[16];
  out.Plant_given = input[17];
  out.Fungal_given = input[18];
  return(out);
};
