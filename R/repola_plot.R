repola_plot <- function() {

  ###
  # Diameter
  ###

  out_needle_mass <- out_m_N_tot <- out_m_R_tot <- out_m_N <- NULL
  count = 1
  n_m_range = seq(0.01, 1, by = 0.1)
  for (n_m in n_m_range) {
    parameters_p_changing <- parameters_p
    parameters_p_changing[c("D0"),1] = n_m

    out_needle_mass[count] = repola_test_cpp(t(parameters_p_changing),
                                              t(sperling_p),
                                              sperling_extras)$needle_mass
    out_m_N_tot[count] = repola_test_cpp(t(parameters_p_changing),
                                          t(sperling_p),
                                          sperling_extras)$m_N_tot
    out_m_R_tot[count] = repola_test_cpp(t(parameters_p_changing),
                                          t(sperling_p),
                                          sperling_extras)$m_R_tot
    out_m_N[count] = repola_test_cpp(t(parameters_p_changing),
                                      t(sperling_p),
                                      sperling_extras)$m_N
    count = count + 1
  }

  plot(n_m_range, out_needle_mass, type = "l", col = "blue",
       main = "Repola Outputs", xlab = "Diameter mass", ylab = "Needle Mass", # TODO: Should it be diameter mass?
       ylim = c(0, max(out_needle_mass)))
  lines(n_m_range, out_m_N_tot, col = "black")
  lines(n_m_range, out_m_R_tot, col = "green")
  lines(n_m_range, out_m_N, col = "purple")

  # TODO: which parameters should be plotted here?
  # TODO: check the paper?
}
