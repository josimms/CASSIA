bound_checks <- function(CASSIA_sensi) {
  direct <- "~/Documents/CASSIA_Calibration/"
  bounds_all <- read.delim(paste0(direct, "bounds_all.csv"), sep = ",", row.names = 1)
  bounds_all$UL[1:20] = rep(10, 20)


  CASSIA_sensi = CASSIA_sensitivity(bounds_all[1:27,],
                                    rownames(bounds_all)[1:27],
                                    2015, 2017, weather_original, GPP_ref,
                                    c(pPREL, N_parameters), t(parameters_p), common_p, t(ratios_p), t(sperling_p),
                                    needle_mass_in,
                                    Throughfall,
                                    storage_rest, storage_grows,
                                    LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                                    trees_grow, growth_decreases, needle_mass_grows,
                                    mycorrhiza, root_as_Ding, sperling_sugar_model,
                                    xylogensis_option, environmental_effect_xylogenesis,
                                    temp_rise, drought, Rm_acclimation,
                                    using_spp_photosynthesis, TRUE,
                                    etmodel, LOGFLAG)

  variables_original <- c("bud", "wall_daily", "needle_daily", "root_daily", "height_daily", "Rg", "Rm", "P") # TODO: check if I want more, and that these are equivalent!
  variables_new <- c("bud_growth", "diameter_growth", "needle_growth", "root_growth", "height_growth", "respiration_growth", "respiration_maintenance", "GPP")

  dates = seq(as.Date("2015-01-01"), as.Date("2017-12-31"), by = "day")
  dates = dates[-c(365+366)]

  Hyde_daily_original <- Hyde_daily_original[1:(365+365+365),]

  ####
  # Growth
  ####

  for (j in 1:27) {
    par(mfrow = c(3, 3), xpd=TRUE)
    for (var in 1:(length(variables_new)+1)) {
      if (var < length(variables_new)) {
        plot(dates, Hyde_daily_original[,c(variables_original[var])],
             main = "", xlab = "Date", ylab = variables_new[var], type = "l")
        for (i in 1:20) {
          lines(dates, CASSIA_sensi[[j]][[i]][[1]][,c(variables_new[var])], col = rainbow(25)[i])
        }
        if (var == 1) {
          title(main = rownames(bounds_all)[j], outer = T, line = -2)
        }
      } else if (var == length(variables_new)) {
        plot(dates, weather_original$P,
             main = "", xlab = "Date", ylab = variables_new[var], type = "l")
        for (i in 1:20) {
          lines(dates, CASSIA_sensi[[j]][[i]][[3]][,c(variables_new[var])], col = rainbow(25)[i])
        }
      } else {
        plot.new()
        legend("center", legend = CASSIA_sensi[[j]][[i]][[4]], col = rainbow(25), lty = 1, bty = "n",
               title = rownames(bounds_all)[j], ncol=2)
      }
    }

    ####
    # Sugar
    ####

    par(mfrow = c(2, 5))
    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$sugar_needles, col = rainbow(25)[1],
         ylab = "Sugar, Needles", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$sugar_needles, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$sugar_phloem, col = rainbow(25)[1],
         ylab = "Sugar, Pholoem", xlab = "Dates", type = "l")
    text(30, 0.43, "Expected Equilibrium", col = "blue", cex = 0.75)
    text(25, sperling_p[c("SCb"),c("Hyde")] + 0.02, "\"bloom\" threshold", col = "pink", cex = 0.75)
    abline(h = 0, lty = 2, col = "grey")
    abline(h = 0.41, lty = 2, col = "blue")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$sugar_phloem, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$sugar_xylem_sh, col = rainbow(25)[1],
         ylab = "Sugar, Xylem sh", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$sugar_xylem_sh, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$sugar_xylem_st, col = rainbow(25)[1],
         ylab = "Sugar, Xylem st", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$sugar_xylem_st, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[i]]$Sugar$sugar_roots, col = rainbow(25)[1],
         ylab = "Sugar, Roots", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$sugar_roots, col = rainbow(25)[i])
    }

    ####
    # Starch
    ####

    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$starch_needles, col = rainbow(25)[1],
         ylab = "Starch, Needles", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$starch_needles, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$starch_phloem, col = rainbow(25)[1],
         ylab = "Starch, Pholoem", xlab = "Dates", type = "l")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$starch_phloem, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$starch_xylem_sh, col = rainbow(25)[1],
         ylab = "Starch, Xylem sh", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$starch_xylem_sh, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[1]]$Sugar$starch_xylem_st, col = rainbow(25)[1],
         ylab = "Starch, Xylem st", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$starch_xylem_st, col = rainbow(25)[i])
    }

    plot(dates, CASSIA_sensi[[j]][[i]]$Sugar$starch_roots, col = rainbow(25)[1],
         ylab = "Starch, Roots", xlab = "Dates", type = "l")
    abline(h = 0, lty = 2, col = "grey")
    for (i in 2:20) {
      lines(dates, CASSIA_sensi[[j]][[i]]$Sugar$starch_roots, col = rainbow(25)[i])
    }
    title(main = rownames(bounds_all)[j], outer = T, line = -2)
  }
}


dhtn <- function(x, mean = 0, sd = 1, log = T) {
  diff <- x-mean
  diff[which(abs(diff)<=1e-9)] <- 1e-9 # to avoid inf when dif <- 0
  R2 <- (diff/sd)^2
  prob <- 1/(sd*(pi*2)^0.5)*(1-exp(-R2/2))/R2
  if (log==F) {return(prob)}
  else {return(log(prob))}
}

likelyhood_sugar_model <- function(par, sum = T) {

  ### Import data
  direct <- "~/Documents/CASSIA_Calibration/"
  load(paste0(direct, "original.data.RData"))
  # All of the observations made into a vector
  # Data processing done outside of the likelihood function to make it quicker
  load(paste0(direct, "obs.vec.og.RData"))

  ### Initial conditions for the CASSIA state values
  # TODO: initial conditions correct!
  sperling_par <- sperling_p
  sperling_par[c(50:54, 35:39, 40:44, 45:49, 25:26),1] <- par[1:22]
  parameters_par <- parameters_p
  parameters_test[c("lower_bound_needles", "lower_bound_phloem", "lower_bound_roots", "lower_bound_xylem_sh", "lower_bound_xylem_st"),1] <- c(0.05, 0.13, 0.007, 0.009, 0.001)
  parameters_par[62:66,1] <- par[23:27]

  sperling_sugar_model = T
  using_spp_photosynthesis = T

  ## TODO: doesn't the function need these values

  simTab <- CASSIA_yearly(2015, 2016, weather_original, GPP_ref,
                          c(pPREL, N_parameters), t(parameters_p), common_p, t(ratios_p), t(sperling_par),
                          needle_mass_in,
                          Throughfall,
                          storage_rest, storage_grows,
                          LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                          trees_grow, growth_decreases, needle_mass_grows,
                          mycorrhiza, root_as_Ding, sperling_sugar_model,
                          xylogensis_option, environmental_effect_xylogenesis,
                          temp_rise, drought, Rm_acclimation,
                          using_spp_photosynthesis, TRUE,
                          etmodel, LOGFLAG)
  simTab_growth = simTab[[1]]
  simTab_sugar = simTab[[2]]
  Hyde_daily_original_cali <- Hyde_daily_original[c(1:730),]

  rownames(simTab_growth) <- rownames(simTab_sugar) <- seq(as.Date("2015-01-01"), as.Date("2016-12-31"), by = "day")[-c(365+366)]
  simTab_sugar <- simTab_sugar[rownames(original.data),]
  # simTab has the same rows as the calibration data, so can directly use the indices
  simVec <- c(simTab_sugar$starch_needles[is.na(original.data$needles_starch) == F],
              simTab_sugar$starch_phloem[is.na(original.data$phloem_starch) == F],
              simTab_sugar$starch_xylem_sh[is.na(original.data$xylem_sh_starch) == F],
              simTab_sugar$starch_xylem_st[is.na(original.data$xylem_st_starch) == F],
              simTab_sugar$starch_roots[is.na(original.data$roots_starch) == F],
              simTab_sugar$sugar_needles[is.na(original.data$needles_sugar) == F],
              simTab_sugar$sugar_phloem[is.na(original.data$phloem_sugar) == F],
              simTab_sugar$sugar_xylem_sh[is.na(original.data$xylem_sh_sugar) == F],
              simTab_sugar$sugar_xylem_st[is.na(original.data$xylem_st_sugar) == F],
              simTab_sugar$sugar_roots[is.na(original.data$roots_sugar) == F],
              simTab_growth$bud_growth,
              simTab_growth$diameter_growth,
              simTab_growth$height_growth,
              simTab_growth$needle_growth,
              simTab_growth$root_growth)
  simVec[simVec < 0 | is.na(simVec) | is.infinite(simVec)] <- 1e6


  ### Standard deviation formulas
  # TODO: look at the residuals of the growth data when it is replaced with the biomass data
  par = bounds_all[,2]
  sdX <- c(par[28]*simVec[1:14]+par[29],
           par[30]*simVec[14+1:14]+par[31],
           par[32]*simVec[28+1:14]+par[33],
           par[34]*simVec[42+1:14]+par[35],
           par[36]*simVec[56+1:13]+par[37],
           par[38]*simVec[69+1:14]+par[39],
           par[40]*simVec[83+1:14]+par[41],
           par[42]*simVec[97+1:14]+par[43],
           par[44]*simVec[111+1:14]+par[45], # TODO Some of the values here have just been set to one - need to decide the variation here!
           par[46]*simVec[125++1:13]+par[47],
           par[48]*simVec[138+1:730]+par[49],
           par[50]*simVec[868+1:730]+par[51],
           par[52]*simVec[1598+1:730]+par[53],
           par[54]*simVec[2328+1:730]+par[55],
           par[56]*simVec[3058+1:730]+par[57])

  res <- c(obs.vec.og[-c(70, 140)],
           Hyde_daily_original_cali$bud,
           Hyde_daily_original_cali$wall_daily,
           Hyde_daily_original_cali$height_daily,
           Hyde_daily_original_cali$needle_daily,
           Hyde_daily_original_cali$root_daily) - simVec

  likeX <- sum(dnorm(res,mean=0,sd=sdX,log=T)) # the likelihood of the model given this parameter

  return(likeX)
}

CASSIA_calibration <- function(preform_callibration = FALSE) {

  # TODO: make new bounds!
  # Run with different parameters sets and check the likelihood, if changes can then use the Bayesian methods
  direct <- "~/Documents/CASSIA_Calibration/"
  bounds_all <- read.delim(paste0(direct, "bounds_all.csv"), sep = ",", row.names = 1)
  bounds_all$UL[1:27] = c(4, 14, 12, 10, 2,
                          12, 12, 12, 12, 12,
                          12, 12, 12, 12, 12,
                          12, 12, 12, 12, 12,
                          6, 6,
                          6, 6, 6, 6, 10)

  # create priors
  prior <- BayesianTools::createUniformPrior(c(bounds_all[,1]),
                                             c(bounds_all[,2]),
                                             best = NULL)

  # Bayesian set up
  CASSIABayesianSetup <- BayesianTools::createBayesianSetup(likelihood = likelyhood_sugar_model,
                                                            prior = prior,
                                                            names = rownames(bounds_all),
                                                            parallel = F)

  settings = list(iterations = 1e6, thin = 1000, nrChains = 3, message = T)

  # Running the MCMC alogorithm
  # Use DEzs as I have read the paper for this one Past samples version
  # If a starting value is given rows are the number of chains and columns are the parameters
  system.time({CASSIAout_sugar_model <- BayesianTools::runMCMC(bayesianSetup = CASSIABayesianSetup, sampler = "DREAMzs", settings = settings)})

  save(CASSIAout_sugar_model, file = paste0(direct, gsub(":", "_", Sys.time()), " CASSIAout_sugar_model.RData"))

  plot(CASSIAout_sugar_model)

  Calibrated_Parameters = BayesianTools::MAP(CASSIAout_sugar_model)$parametersMAP

  test_against_original_data(c(Calibrated_Parameters), T)

  return(CASSIAout_sugar_model)
}

