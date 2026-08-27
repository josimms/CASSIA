#include "CASSIA.h"


int CASSIA_ratios_test(std::string fratios, std::string site) {
  // Example usage: read values for the 'Hyde' site
  CASSIA_ratios hyde_ratios = read_ratios(fratios, site);

  // Output the read values
  std::cout << "Hyde Ratios:" << std::endl;
  std::cout << "Form Factor: " << hyde_ratios.form_factor << std::endl;
  std::cout << "Needle Fine Root Ratio: " << hyde_ratios.needle_fineroot_ratio << std::endl;
  std::cout << "Sapwood Share: " << hyde_ratios.sapwood_share << std::endl;
  std::cout << "Height Growth Coefficient: " << hyde_ratios.height_growth_coefficient << std::endl;
  std::cout << "Diameter Growth Coefficient: " << hyde_ratios.diameter_growth_coefficient << std::endl;
  std::cout << "Height Growth Coefficient Max: " << hyde_ratios.height_growth_coefficient_max << std::endl;
  std::cout << "Height Growth Coefficient Min: " << hyde_ratios.height_growth_coefficient_min << std::endl;
  std::cout << "Diameter Growth Coefficient Max: " << hyde_ratios.diameter_growth_coefficient_max << std::endl;
  std::cout << "Diameter Growth Coefficient Min: " << hyde_ratios.diameter_growth_coefficient_min << std::endl;

  return 0;
}


int CASSIA_common_test(std::string fcommon) {
  CASSIA_common common_params = read_common_parameters(fcommon);

  // Output the read values
  std::cout << "CASSIA Common Parameters:" << std::endl;
  std::cout << "a: " << std::setprecision(16) << common_params.a << std::endl;
  std::cout << "b: " << common_params.b << std::endl;
  std::cout << "TR0: " << common_params.TR0 << std::endl;
  std::cout << "abs_zero: " << common_params.abs_zero << std::endl;
  std::cout << "b_s: " << common_params.b_s << std::endl;
  std::cout << "theetta_FC: " << common_params.theetta_FC << std::endl;
  std::cout << "phi_e: " << common_params.phi_e << std::endl;
  std::cout << "K_sat: " << common_params.K_sat << std::endl;
  std::cout << "R_length: " << common_params.R_length << std::endl;
  std::cout << "M_H20: " << common_params.M_H20 << std::endl;
  std::cout << "r_cyl: " << common_params.r_cyl << std::endl;
  std::cout << "r_root: " << common_params.r_root << std::endl;
  std::cout << "ypsilon: " << common_params.ypsilon << std::endl;
  std::cout << "Rg_N: " << common_params.Rg_N << std::endl;
  std::cout << "Rg_S: " << common_params.Rg_S << std::endl;
  std::cout << "Rg_R: " << common_params.Rg_R << std::endl;
  std::cout << "gas_const: " << common_params.gas_const << std::endl;
  std::cout << "M_C: " << common_params.M_C << std::endl;
  std::cout << "M_H: " << common_params.M_H << std::endl;
  std::cout << "M_O: " << common_params.M_O << std::endl;
  std::cout << "osmotic_sugar_conc: " << common_params.osmotic_sugar_conc << std::endl;
  std::cout << "m_N: " << common_params.m_N << std::endl;
  std::cout << "Uggla: " << common_params.Uggla << std::endl;

  return 0;
}


int CASSIA_parameter_test(CASSIA_parameters structure_params) {

#define PRINT_PARAM(param) std::cout << #param << ": " << structure_params.param << std::endl

  PRINT_PARAM(Q10_N);
  PRINT_PARAM(Rm_N);
  PRINT_PARAM(Q10_S);
  PRINT_PARAM(Rm_S);
  PRINT_PARAM(Q10_R);
  PRINT_PARAM(Rm_R);
  PRINT_PARAM(sR0);
  PRINT_PARAM(sRc);
  PRINT_PARAM(growth_myco);
  PRINT_PARAM(root_lifetime);
  PRINT_PARAM(HH0);
  PRINT_PARAM(sH0);
  PRINT_PARAM(LH);
  PRINT_PARAM(LH0);
  PRINT_PARAM(sHc);
  PRINT_PARAM(sN0);
  PRINT_PARAM(LN);
  PRINT_PARAM(LN0);
  PRINT_PARAM(sNc);
  PRINT_PARAM(HN0);
  PRINT_PARAM(sD0_Trad);
  PRINT_PARAM(LD);
  PRINT_PARAM(LD0);
  PRINT_PARAM(sDc);
  PRINT_PARAM(sDc_T_count);
  PRINT_PARAM(tau_Ee);
  PRINT_PARAM(tau_El);
  PRINT_PARAM(tau_We);
  PRINT_PARAM(tau_Wl);
  PRINT_PARAM(tau_GPP);
  PRINT_PARAM(Uggla);
  PRINT_PARAM(sB0);
  PRINT_PARAM(sBc);
  PRINT_PARAM(LB);
  PRINT_PARAM(cell_d_ew);
  PRINT_PARAM(cell_d_lw);
  PRINT_PARAM(cell_l_ew);
  PRINT_PARAM(cell_l_lw);
  PRINT_PARAM(cell_wall_density_ew);
  PRINT_PARAM(cell_wall_density_lw);
  PRINT_PARAM(wall_thickness_ew);
  PRINT_PARAM(wall_thickness_lw);
  PRINT_PARAM(cell_volume_growth_per_day_ew);
  PRINT_PARAM(cell_volume_growth_per_day_lw);
  PRINT_PARAM(density_tree);
  PRINT_PARAM(carbon_share);
  PRINT_PARAM(D0);
  PRINT_PARAM(h0);
  PRINT_PARAM(n_age);
  PRINT_PARAM(n_length);
  PRINT_PARAM(h_increment);
  PRINT_PARAM(SLA);
  PRINT_PARAM(LR0);
  PRINT_PARAM(mycorrhiza_threshold);
  PRINT_PARAM(starch0);
  PRINT_PARAM(sugar0);
  PRINT_PARAM(starch_needles0);
  PRINT_PARAM(starch_phloem0);
  PRINT_PARAM(starch_xylem_sh0);
  PRINT_PARAM(starch_xylem_st0);
  PRINT_PARAM(starch_roots0);
  PRINT_PARAM(sugar_needles0);
  PRINT_PARAM(sugar_phloem0);
  PRINT_PARAM(sugar_roots0);
  PRINT_PARAM(sugar_xylem_sh0);
  PRINT_PARAM(sugar_xylem_st0);
  PRINT_PARAM(Wala_needles);
  PRINT_PARAM(Wala_phloem);
  PRINT_PARAM(Wala_xylem_sh);
  PRINT_PARAM(Wala_xylem_st);
  PRINT_PARAM(Wala_roots);
  PRINT_PARAM(carbon_sugar);
  PRINT_PARAM(carbon_starch);
  PRINT_PARAM(alfa_needles);
  PRINT_PARAM(alfa_phloem);
  PRINT_PARAM(alfa_xylem_sh);
  PRINT_PARAM(alfa_xylem_st);
  PRINT_PARAM(alfa_roots);
  PRINT_PARAM(tau_s);
  PRINT_PARAM(tau_t);
  PRINT_PARAM(starch00);
  PRINT_PARAM(sugar00);
  PRINT_PARAM(Wala);
  PRINT_PARAM(alfa);
  PRINT_PARAM(Q10s);
  PRINT_PARAM(Q10d);
  PRINT_PARAM(SCb);
  PRINT_PARAM(sugar_level);
  PRINT_PARAM(Ad0_needles);
  PRINT_PARAM(Ad0_phloem);
  PRINT_PARAM(Ad0_roots);
  PRINT_PARAM(Ad0_xylem_sh);
  PRINT_PARAM(Ad0_xylem_st);
  PRINT_PARAM(lambda_needles);
  PRINT_PARAM(lambda_phloem);
  PRINT_PARAM(lambda_roots);
  PRINT_PARAM(lambda_xylem_sh);
  PRINT_PARAM(lambda_xylem_st);
  PRINT_PARAM(delta_needles);
  PRINT_PARAM(delta_phloem);
  PRINT_PARAM(delta_roots);
  PRINT_PARAM(delta_xylem_sh);
  PRINT_PARAM(delta_xylem_st);
  PRINT_PARAM(lower_bound_needles);
  PRINT_PARAM(lower_bound_phloem);
  PRINT_PARAM(lower_bound_roots);
  PRINT_PARAM(lower_bound_xylem_sh);
  PRINT_PARAM(lower_bound_xylem_st);
  PRINT_PARAM(tau_emergancy_needles);
  PRINT_PARAM(tau_emergancy_phloem);
  PRINT_PARAM(tau_emergancy_roots);
  PRINT_PARAM(tau_emergancy_xylem_sh);
  PRINT_PARAM(tau_emergancy_xylem_st);
  PRINT_PARAM(percentage_needle_storage);
  PRINT_PARAM(percentage_phloem_storage);
  PRINT_PARAM(percentage_xylem_sh_storage);
  PRINT_PARAM(percentage_xylem_st_storage);
  PRINT_PARAM(percentage_roots_storage);
  PRINT_PARAM(lower_bound_W);
  PRINT_PARAM(tau_emergancy);
  PRINT_PARAM(b0_repo);
  PRINT_PARAM(b1_repo);
  PRINT_PARAM(b2_repo);
  PRINT_PARAM(uk_repo);
  PRINT_PARAM(eki_repo);
  PRINT_PARAM(stem_no);
  PRINT_PARAM(m_R_tot);
  PRINT_PARAM(diameter_start_day);

  return 0;
}
