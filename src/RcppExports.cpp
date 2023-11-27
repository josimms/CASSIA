// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "CASSIA_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CASSIA_yearly
Rcpp::List CASSIA_yearly(int start_year, int end_year, Rcpp::DataFrame weather, std::vector<double> GPP_ref, std::vector<double> pPREL, Rcpp::DataFrame pCASSIA_parameters, Rcpp::DataFrame pCASSIA_common, Rcpp::DataFrame pCASSIA_ratios, Rcpp::DataFrame pCASSIA_sperling, double needle_mass_in, double Throughfall, bool storage_rest, bool storage_grows, bool LH_estim, bool LN_estim, bool mN_varies, bool LD_estim, bool sD_estim_T_count, bool trees_grow, bool growth_decreases, bool needle_mass_grows, bool mycorrhiza, bool root_as_Ding, bool sperling_sugar_model, bool xylogensis_option, bool environmental_effect_xylogenesis, bool temp_rise, bool drought, bool Rm_acclimation, bool using_spp_photosynthesis, bool CASSIA_graphs, int etmodel, int LOGFLAG);
RcppExport SEXP _CASSIA_CASSIA_yearly(SEXP start_yearSEXP, SEXP end_yearSEXP, SEXP weatherSEXP, SEXP GPP_refSEXP, SEXP pPRELSEXP, SEXP pCASSIA_parametersSEXP, SEXP pCASSIA_commonSEXP, SEXP pCASSIA_ratiosSEXP, SEXP pCASSIA_sperlingSEXP, SEXP needle_mass_inSEXP, SEXP ThroughfallSEXP, SEXP storage_restSEXP, SEXP storage_growsSEXP, SEXP LH_estimSEXP, SEXP LN_estimSEXP, SEXP mN_variesSEXP, SEXP LD_estimSEXP, SEXP sD_estim_T_countSEXP, SEXP trees_growSEXP, SEXP growth_decreasesSEXP, SEXP needle_mass_growsSEXP, SEXP mycorrhizaSEXP, SEXP root_as_DingSEXP, SEXP sperling_sugar_modelSEXP, SEXP xylogensis_optionSEXP, SEXP environmental_effect_xylogenesisSEXP, SEXP temp_riseSEXP, SEXP droughtSEXP, SEXP Rm_acclimationSEXP, SEXP using_spp_photosynthesisSEXP, SEXP CASSIA_graphsSEXP, SEXP etmodelSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type start_year(start_yearSEXP);
    Rcpp::traits::input_parameter< int >::type end_year(end_yearSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type weather(weatherSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type GPP_ref(GPP_refSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pPREL(pPRELSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_parameters(pCASSIA_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_common(pCASSIA_commonSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_ratios(pCASSIA_ratiosSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_sperling(pCASSIA_sperlingSEXP);
    Rcpp::traits::input_parameter< double >::type needle_mass_in(needle_mass_inSEXP);
    Rcpp::traits::input_parameter< double >::type Throughfall(ThroughfallSEXP);
    Rcpp::traits::input_parameter< bool >::type storage_rest(storage_restSEXP);
    Rcpp::traits::input_parameter< bool >::type storage_grows(storage_growsSEXP);
    Rcpp::traits::input_parameter< bool >::type LH_estim(LH_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type LN_estim(LN_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type mN_varies(mN_variesSEXP);
    Rcpp::traits::input_parameter< bool >::type LD_estim(LD_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type sD_estim_T_count(sD_estim_T_countSEXP);
    Rcpp::traits::input_parameter< bool >::type trees_grow(trees_growSEXP);
    Rcpp::traits::input_parameter< bool >::type growth_decreases(growth_decreasesSEXP);
    Rcpp::traits::input_parameter< bool >::type needle_mass_grows(needle_mass_growsSEXP);
    Rcpp::traits::input_parameter< bool >::type mycorrhiza(mycorrhizaSEXP);
    Rcpp::traits::input_parameter< bool >::type root_as_Ding(root_as_DingSEXP);
    Rcpp::traits::input_parameter< bool >::type sperling_sugar_model(sperling_sugar_modelSEXP);
    Rcpp::traits::input_parameter< bool >::type xylogensis_option(xylogensis_optionSEXP);
    Rcpp::traits::input_parameter< bool >::type environmental_effect_xylogenesis(environmental_effect_xylogenesisSEXP);
    Rcpp::traits::input_parameter< bool >::type temp_rise(temp_riseSEXP);
    Rcpp::traits::input_parameter< bool >::type drought(droughtSEXP);
    Rcpp::traits::input_parameter< bool >::type Rm_acclimation(Rm_acclimationSEXP);
    Rcpp::traits::input_parameter< bool >::type using_spp_photosynthesis(using_spp_photosynthesisSEXP);
    Rcpp::traits::input_parameter< bool >::type CASSIA_graphs(CASSIA_graphsSEXP);
    Rcpp::traits::input_parameter< int >::type etmodel(etmodelSEXP);
    Rcpp::traits::input_parameter< int >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(CASSIA_yearly(start_year, end_year, weather, GPP_ref, pPREL, pCASSIA_parameters, pCASSIA_common, pCASSIA_ratios, pCASSIA_sperling, needle_mass_in, Throughfall, storage_rest, storage_grows, LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count, trees_grow, growth_decreases, needle_mass_grows, mycorrhiza, root_as_Ding, sperling_sugar_model, xylogensis_option, environmental_effect_xylogenesis, temp_rise, drought, Rm_acclimation, using_spp_photosynthesis, CASSIA_graphs, etmodel, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// replace_value_DataFrame
Rcpp::DataFrame replace_value_DataFrame(Rcpp::DataFrame df, double value, int ref);
RcppExport SEXP _CASSIA_replace_value_DataFrame(SEXP dfSEXP, SEXP valueSEXP, SEXP refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type ref(refSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_value_DataFrame(df, value, ref));
    return rcpp_result_gen;
END_RCPP
}
// CASSIA_sensitivity
Rcpp::List CASSIA_sensitivity(Rcpp::DataFrame bounds, std::vector<std::string> names, int start_year, int end_year, Rcpp::DataFrame weather, std::vector<double> GPP_ref, std::vector<double> pPREL, Rcpp::DataFrame pCASSIA_parameters, Rcpp::DataFrame pCASSIA_common, Rcpp::DataFrame pCASSIA_ratios, Rcpp::DataFrame pCASSIA_sperling, double needle_mass_in, double Throughfall, bool storage_rest, bool storage_grows, bool LH_estim, bool LN_estim, bool mN_varies, bool LD_estim, bool sD_estim_T_count, bool trees_grow, bool growth_decreases, bool needle_mass_grows, bool mycorrhiza, bool root_as_Ding, bool sperling_sugar_model, bool xylogensis_option, bool environmental_effect_xylogenesis, bool temp_rise, bool drought, bool Rm_acclimation, bool using_spp_photosynthesis, bool CASSIA_graphs, int etmodel, int LOGFLAG);
RcppExport SEXP _CASSIA_CASSIA_sensitivity(SEXP boundsSEXP, SEXP namesSEXP, SEXP start_yearSEXP, SEXP end_yearSEXP, SEXP weatherSEXP, SEXP GPP_refSEXP, SEXP pPRELSEXP, SEXP pCASSIA_parametersSEXP, SEXP pCASSIA_commonSEXP, SEXP pCASSIA_ratiosSEXP, SEXP pCASSIA_sperlingSEXP, SEXP needle_mass_inSEXP, SEXP ThroughfallSEXP, SEXP storage_restSEXP, SEXP storage_growsSEXP, SEXP LH_estimSEXP, SEXP LN_estimSEXP, SEXP mN_variesSEXP, SEXP LD_estimSEXP, SEXP sD_estim_T_countSEXP, SEXP trees_growSEXP, SEXP growth_decreasesSEXP, SEXP needle_mass_growsSEXP, SEXP mycorrhizaSEXP, SEXP root_as_DingSEXP, SEXP sperling_sugar_modelSEXP, SEXP xylogensis_optionSEXP, SEXP environmental_effect_xylogenesisSEXP, SEXP temp_riseSEXP, SEXP droughtSEXP, SEXP Rm_acclimationSEXP, SEXP using_spp_photosynthesisSEXP, SEXP CASSIA_graphsSEXP, SEXP etmodelSEXP, SEXP LOGFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type names(namesSEXP);
    Rcpp::traits::input_parameter< int >::type start_year(start_yearSEXP);
    Rcpp::traits::input_parameter< int >::type end_year(end_yearSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type weather(weatherSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type GPP_ref(GPP_refSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pPREL(pPRELSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_parameters(pCASSIA_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_common(pCASSIA_commonSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_ratios(pCASSIA_ratiosSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_sperling(pCASSIA_sperlingSEXP);
    Rcpp::traits::input_parameter< double >::type needle_mass_in(needle_mass_inSEXP);
    Rcpp::traits::input_parameter< double >::type Throughfall(ThroughfallSEXP);
    Rcpp::traits::input_parameter< bool >::type storage_rest(storage_restSEXP);
    Rcpp::traits::input_parameter< bool >::type storage_grows(storage_growsSEXP);
    Rcpp::traits::input_parameter< bool >::type LH_estim(LH_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type LN_estim(LN_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type mN_varies(mN_variesSEXP);
    Rcpp::traits::input_parameter< bool >::type LD_estim(LD_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type sD_estim_T_count(sD_estim_T_countSEXP);
    Rcpp::traits::input_parameter< bool >::type trees_grow(trees_growSEXP);
    Rcpp::traits::input_parameter< bool >::type growth_decreases(growth_decreasesSEXP);
    Rcpp::traits::input_parameter< bool >::type needle_mass_grows(needle_mass_growsSEXP);
    Rcpp::traits::input_parameter< bool >::type mycorrhiza(mycorrhizaSEXP);
    Rcpp::traits::input_parameter< bool >::type root_as_Ding(root_as_DingSEXP);
    Rcpp::traits::input_parameter< bool >::type sperling_sugar_model(sperling_sugar_modelSEXP);
    Rcpp::traits::input_parameter< bool >::type xylogensis_option(xylogensis_optionSEXP);
    Rcpp::traits::input_parameter< bool >::type environmental_effect_xylogenesis(environmental_effect_xylogenesisSEXP);
    Rcpp::traits::input_parameter< bool >::type temp_rise(temp_riseSEXP);
    Rcpp::traits::input_parameter< bool >::type drought(droughtSEXP);
    Rcpp::traits::input_parameter< bool >::type Rm_acclimation(Rm_acclimationSEXP);
    Rcpp::traits::input_parameter< bool >::type using_spp_photosynthesis(using_spp_photosynthesisSEXP);
    Rcpp::traits::input_parameter< bool >::type CASSIA_graphs(CASSIA_graphsSEXP);
    Rcpp::traits::input_parameter< int >::type etmodel(etmodelSEXP);
    Rcpp::traits::input_parameter< int >::type LOGFLAG(LOGFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(CASSIA_sensitivity(bounds, names, start_year, end_year, weather, GPP_ref, pPREL, pCASSIA_parameters, pCASSIA_common, pCASSIA_ratios, pCASSIA_sperling, needle_mass_in, Throughfall, storage_rest, storage_grows, LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count, trees_grow, growth_decreases, needle_mass_grows, mycorrhiza, root_as_Ding, sperling_sugar_model, xylogensis_option, environmental_effect_xylogenesis, temp_rise, drought, Rm_acclimation, using_spp_photosynthesis, CASSIA_graphs, etmodel, LOGFLAG));
    return rcpp_result_gen;
END_RCPP
}
// xylogenesis_wrapper
Rcpp::List xylogenesis_wrapper(int no_day, int day, Rcpp::DataFrame pCASSIA_parameters, Rcpp::DataFrame pCASSIA_common, Rcpp::DataFrame pCASSIA_sperling, std::vector<double> extras_sperling, bool xylogenesis_option, bool environmental_effect_xylogenesis, double TAir, double n_rows, double max_ew_cells, double n_E_pot_old, double n_W_pot_old, double n_M_pot_old, double g, std::vector<double> en_growth_vector, double tau_W_old, double carbon_daily_rate_ew, double carbon_daily_rate_lw);
RcppExport SEXP _CASSIA_xylogenesis_wrapper(SEXP no_daySEXP, SEXP daySEXP, SEXP pCASSIA_parametersSEXP, SEXP pCASSIA_commonSEXP, SEXP pCASSIA_sperlingSEXP, SEXP extras_sperlingSEXP, SEXP xylogenesis_optionSEXP, SEXP environmental_effect_xylogenesisSEXP, SEXP TAirSEXP, SEXP n_rowsSEXP, SEXP max_ew_cellsSEXP, SEXP n_E_pot_oldSEXP, SEXP n_W_pot_oldSEXP, SEXP n_M_pot_oldSEXP, SEXP gSEXP, SEXP en_growth_vectorSEXP, SEXP tau_W_oldSEXP, SEXP carbon_daily_rate_ewSEXP, SEXP carbon_daily_rate_lwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type no_day(no_daySEXP);
    Rcpp::traits::input_parameter< int >::type day(daySEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_parameters(pCASSIA_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_common(pCASSIA_commonSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_sperling(pCASSIA_sperlingSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type extras_sperling(extras_sperlingSEXP);
    Rcpp::traits::input_parameter< bool >::type xylogenesis_option(xylogenesis_optionSEXP);
    Rcpp::traits::input_parameter< bool >::type environmental_effect_xylogenesis(environmental_effect_xylogenesisSEXP);
    Rcpp::traits::input_parameter< double >::type TAir(TAirSEXP);
    Rcpp::traits::input_parameter< double >::type n_rows(n_rowsSEXP);
    Rcpp::traits::input_parameter< double >::type max_ew_cells(max_ew_cellsSEXP);
    Rcpp::traits::input_parameter< double >::type n_E_pot_old(n_E_pot_oldSEXP);
    Rcpp::traits::input_parameter< double >::type n_W_pot_old(n_W_pot_oldSEXP);
    Rcpp::traits::input_parameter< double >::type n_M_pot_old(n_M_pot_oldSEXP);
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type en_growth_vector(en_growth_vectorSEXP);
    Rcpp::traits::input_parameter< double >::type tau_W_old(tau_W_oldSEXP);
    Rcpp::traits::input_parameter< double >::type carbon_daily_rate_ew(carbon_daily_rate_ewSEXP);
    Rcpp::traits::input_parameter< double >::type carbon_daily_rate_lw(carbon_daily_rate_lwSEXP);
    rcpp_result_gen = Rcpp::wrap(xylogenesis_wrapper(no_day, day, pCASSIA_parameters, pCASSIA_common, pCASSIA_sperling, extras_sperling, xylogenesis_option, environmental_effect_xylogenesis, TAir, n_rows, max_ew_cells, n_E_pot_old, n_W_pot_old, n_M_pot_old, g, en_growth_vector, tau_W_old, carbon_daily_rate_ew, carbon_daily_rate_lw));
    return rcpp_result_gen;
END_RCPP
}
// growth_wrapper
Rcpp::List growth_wrapper(int day, int year, double TAir, double TSoil_A, double TSoil_B, double Soil_Moisture, double PF, double GPP_ref, bool root_as_Ding, bool xylogenesis_option, bool environmental_effect_xylogenesis, bool sD_estim_T_count, Rcpp::DataFrame pCASSIA_common, Rcpp::DataFrame pCASSIA_parameters, Rcpp::DataFrame pCASSIA_ratios, Rcpp::DataFrame pCASSIA_sperling, std::vector<double> extras_sperling, double CH, double B0, double en_pot_growth_old, double GPP_mean, double GPP_previous_sum, bool LH_estim, bool LN_estim, bool LD_estim, std::vector<double> growth_in, double last_year_HH, int no_day);
RcppExport SEXP _CASSIA_growth_wrapper(SEXP daySEXP, SEXP yearSEXP, SEXP TAirSEXP, SEXP TSoil_ASEXP, SEXP TSoil_BSEXP, SEXP Soil_MoistureSEXP, SEXP PFSEXP, SEXP GPP_refSEXP, SEXP root_as_DingSEXP, SEXP xylogenesis_optionSEXP, SEXP environmental_effect_xylogenesisSEXP, SEXP sD_estim_T_countSEXP, SEXP pCASSIA_commonSEXP, SEXP pCASSIA_parametersSEXP, SEXP pCASSIA_ratiosSEXP, SEXP pCASSIA_sperlingSEXP, SEXP extras_sperlingSEXP, SEXP CHSEXP, SEXP B0SEXP, SEXP en_pot_growth_oldSEXP, SEXP GPP_meanSEXP, SEXP GPP_previous_sumSEXP, SEXP LH_estimSEXP, SEXP LN_estimSEXP, SEXP LD_estimSEXP, SEXP growth_inSEXP, SEXP last_year_HHSEXP, SEXP no_daySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type day(daySEXP);
    Rcpp::traits::input_parameter< int >::type year(yearSEXP);
    Rcpp::traits::input_parameter< double >::type TAir(TAirSEXP);
    Rcpp::traits::input_parameter< double >::type TSoil_A(TSoil_ASEXP);
    Rcpp::traits::input_parameter< double >::type TSoil_B(TSoil_BSEXP);
    Rcpp::traits::input_parameter< double >::type Soil_Moisture(Soil_MoistureSEXP);
    Rcpp::traits::input_parameter< double >::type PF(PFSEXP);
    Rcpp::traits::input_parameter< double >::type GPP_ref(GPP_refSEXP);
    Rcpp::traits::input_parameter< bool >::type root_as_Ding(root_as_DingSEXP);
    Rcpp::traits::input_parameter< bool >::type xylogenesis_option(xylogenesis_optionSEXP);
    Rcpp::traits::input_parameter< bool >::type environmental_effect_xylogenesis(environmental_effect_xylogenesisSEXP);
    Rcpp::traits::input_parameter< bool >::type sD_estim_T_count(sD_estim_T_countSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_common(pCASSIA_commonSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_parameters(pCASSIA_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_ratios(pCASSIA_ratiosSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_sperling(pCASSIA_sperlingSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type extras_sperling(extras_sperlingSEXP);
    Rcpp::traits::input_parameter< double >::type CH(CHSEXP);
    Rcpp::traits::input_parameter< double >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< double >::type en_pot_growth_old(en_pot_growth_oldSEXP);
    Rcpp::traits::input_parameter< double >::type GPP_mean(GPP_meanSEXP);
    Rcpp::traits::input_parameter< double >::type GPP_previous_sum(GPP_previous_sumSEXP);
    Rcpp::traits::input_parameter< bool >::type LH_estim(LH_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type LN_estim(LN_estimSEXP);
    Rcpp::traits::input_parameter< bool >::type LD_estim(LD_estimSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type growth_in(growth_inSEXP);
    Rcpp::traits::input_parameter< double >::type last_year_HH(last_year_HHSEXP);
    Rcpp::traits::input_parameter< int >::type no_day(no_daySEXP);
    rcpp_result_gen = Rcpp::wrap(growth_wrapper(day, year, TAir, TSoil_A, TSoil_B, Soil_Moisture, PF, GPP_ref, root_as_Ding, xylogenesis_option, environmental_effect_xylogenesis, sD_estim_T_count, pCASSIA_common, pCASSIA_parameters, pCASSIA_ratios, pCASSIA_sperling, extras_sperling, CH, B0, en_pot_growth_old, GPP_mean, GPP_previous_sum, LH_estim, LN_estim, LD_estim, growth_in, last_year_HH, no_day));
    return rcpp_result_gen;
END_RCPP
}
// preles_test_cpp
Rcpp::List preles_test_cpp(int NofDays, int day, Rcpp::DataFrame weather, std::vector<double> pPREL, int etmodel);
RcppExport SEXP _CASSIA_preles_test_cpp(SEXP NofDaysSEXP, SEXP daySEXP, SEXP weatherSEXP, SEXP pPRELSEXP, SEXP etmodelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type NofDays(NofDaysSEXP);
    Rcpp::traits::input_parameter< int >::type day(daySEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type weather(weatherSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pPREL(pPRELSEXP);
    Rcpp::traits::input_parameter< int >::type etmodel(etmodelSEXP);
    rcpp_result_gen = Rcpp::wrap(preles_test_cpp(NofDays, day, weather, pPREL, etmodel));
    return rcpp_result_gen;
END_RCPP
}
// repola_test_cpp
Rcpp::List repola_test_cpp(Rcpp::DataFrame pCASSIA_parameters, Rcpp::DataFrame pCASSIA_sperling);
RcppExport SEXP _CASSIA_repola_test_cpp(SEXP pCASSIA_parametersSEXP, SEXP pCASSIA_sperlingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_parameters(pCASSIA_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_sperling(pCASSIA_sperlingSEXP);
    rcpp_result_gen = Rcpp::wrap(repola_test_cpp(pCASSIA_parameters, pCASSIA_sperling));
    return rcpp_result_gen;
END_RCPP
}
// respiration_test_cpp
Rcpp::List respiration_test_cpp(Rcpp::DataFrame pCASSIA_parameters, Rcpp::DataFrame pCASSIA_common, Rcpp::DataFrame pCASSIA_ratios, Rcpp::DataFrame pCASSIA_sperling, std::vector<double> extras_sperling, int day, double TAir, double TSoil, bool temp_rise, bool Rm_acclimation, bool mN_varies, double B0);
RcppExport SEXP _CASSIA_respiration_test_cpp(SEXP pCASSIA_parametersSEXP, SEXP pCASSIA_commonSEXP, SEXP pCASSIA_ratiosSEXP, SEXP pCASSIA_sperlingSEXP, SEXP extras_sperlingSEXP, SEXP daySEXP, SEXP TAirSEXP, SEXP TSoilSEXP, SEXP temp_riseSEXP, SEXP Rm_acclimationSEXP, SEXP mN_variesSEXP, SEXP B0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_parameters(pCASSIA_parametersSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_common(pCASSIA_commonSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_ratios(pCASSIA_ratiosSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type pCASSIA_sperling(pCASSIA_sperlingSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type extras_sperling(extras_sperlingSEXP);
    Rcpp::traits::input_parameter< int >::type day(daySEXP);
    Rcpp::traits::input_parameter< double >::type TAir(TAirSEXP);
    Rcpp::traits::input_parameter< double >::type TSoil(TSoilSEXP);
    Rcpp::traits::input_parameter< bool >::type temp_rise(temp_riseSEXP);
    Rcpp::traits::input_parameter< bool >::type Rm_acclimation(Rm_acclimationSEXP);
    Rcpp::traits::input_parameter< bool >::type mN_varies(mN_variesSEXP);
    Rcpp::traits::input_parameter< double >::type B0(B0SEXP);
    rcpp_result_gen = Rcpp::wrap(respiration_test_cpp(pCASSIA_parameters, pCASSIA_common, pCASSIA_ratios, pCASSIA_sperling, extras_sperling, day, TAir, TSoil, temp_rise, Rm_acclimation, mN_varies, B0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CASSIA_CASSIA_yearly", (DL_FUNC) &_CASSIA_CASSIA_yearly, 33},
    {"_CASSIA_replace_value_DataFrame", (DL_FUNC) &_CASSIA_replace_value_DataFrame, 3},
    {"_CASSIA_CASSIA_sensitivity", (DL_FUNC) &_CASSIA_CASSIA_sensitivity, 35},
    {"_CASSIA_xylogenesis_wrapper", (DL_FUNC) &_CASSIA_xylogenesis_wrapper, 19},
    {"_CASSIA_growth_wrapper", (DL_FUNC) &_CASSIA_growth_wrapper, 28},
    {"_CASSIA_preles_test_cpp", (DL_FUNC) &_CASSIA_preles_test_cpp, 5},
    {"_CASSIA_repola_test_cpp", (DL_FUNC) &_CASSIA_repola_test_cpp, 2},
    {"_CASSIA_respiration_test_cpp", (DL_FUNC) &_CASSIA_respiration_test_cpp, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_CASSIA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}