xylem_plot <- function(pCASSIA_parameters, pCASSIA_common, pCASSIA_sperling, extras_sperling,
                       xylogenesis_option, environmental_effect_xylogenesis, data_format,
                       n_rows, max_ew_cells,
                       n_E_pot_old, n_W_pot_old, n_M_pot_old, g, en_growth_vector,
                       tau_W_old, carbon_daily_rate_ew, carbon_daily_rate_lw) {

  needle_mass <- m_N_tot <- m_R_tot <- m_N <- NULL
  for (day in 1:nrow(data_format)) {
    needle_mass[day] = xylogenesis_wrapper(nrow(data_format), day,
                                           pCASSIA_parameters, pCASSIA_common, pCASSIA_sperling, extras_sperling, # TODO: get the parameters from somewhere
                                           xylogenesis_option, environmental_effect_xylogenesis,
                                           data_format$T[day], n_rows, max_ew_cells,
                                           n_E_pot_old, n_W_pot_old, n_M_pot_old, g, en_growth_vector,
                                           tau_W_old, carbon_daily_rate_ew, carbon_daily_rate_lw)$needle_mass
    m_N_tot[day] = xylogenesis_wrapper(nrow(data_format), day,
                                       pCASSIA_parameters, pCASSIA_common, pCASSIA_sperling, extras_sperling, # TODO: get the parameters from somewhere
                                       xylogenesis_option, environmental_effect_xylogenesis,
                                       data_format$T[day], n_rows, max_ew_cells,
                                       n_E_pot_old, n_W_pot_old, n_M_pot_old, g, en_growth_vector,
                                       tau_W_old, carbon_daily_rate_ew, carbon_daily_rate_lw)$m_N_tot
    m_R_tot[day] = xylogenesis_wrapper(nrow(data_format), day,
                                       pCASSIA_parameters, pCASSIA_common, pCASSIA_sperling, extras_sperling, # TODO: get the parameters from somewhere
                                       xylogenesis_option, environmental_effect_xylogenesis,
                                       data_format$T[day], n_rows, max_ew_cells,
                                       n_E_pot_old, n_W_pot_old, n_M_pot_old, g, en_growth_vector,
                                       tau_W_old, carbon_daily_rate_ew, carbon_daily_rate_lw)$m_R_tot
    m_N[day] = xylogenesis_wrapper(nrow(data_format), day,
                                   pCASSIA_parameters, pCASSIA_common, pCASSIA_sperling, extras_sperling, # TODO: get the parameters from somewhere
                                   xylogenesis_option, environmental_effect_xylogenesis,
                                   data_format$T[day], n_rows, max_ew_cells,
                                   n_E_pot_old, n_W_pot_old, n_M_pot_old, g, en_growth_vector,
                                   tau_W_old, carbon_daily_rate_ew, carbon_daily_rate_lw)$m_N
  }

  par(mfrow = c(2, 2))
  plot(data_format$Date, needle_mass, main = "Needle Mass", xlab = "Date", ylab = "")
  plot(data_format$Date, m_N_tot, main = "m_N_tot", xlab = "Date", ylab = "")
  plot(data_format$Date, m_R_tot, main = "m_R_tot", xlab = "Date", ylab = "")
  plot(data_format$Date, m_N, main = "m_N", xlab = "Date", ylab = "")
}
