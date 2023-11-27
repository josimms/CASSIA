#' @export
sugar_plot <- function(CASSIA_out, site, SCb = sperling_p[c("SCb"), c(site)], sugar_model_cpp_out, cpp = T) {

  if (cpp) {
    par(mfrow = c(2, 1))
    plot(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$sugar_needles +
           sugar_model_cpp$Sugar$sugar_phloem +
           sugar_model_cpp$Sugar$sugar_roots +
           sugar_model_cpp$Sugar$sugar_xylem_sh +
           sugar_model_cpp$Sugar$sugar_xylem_st, main = "Sugar", col = "blue", type = "l", ylim = c(0, 0.8), xlab = "Days of the Year", ylab = "Sugar, kg C")
    abline(h = 0, lty = 2, col = "grey")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$sugar_phloem +
            sugar_model_cpp$Sugar$sugar_roots +
            sugar_model_cpp$Sugar$sugar_xylem_sh +
            sugar_model_cpp$Sugar$sugar_xylem_st, col = "blue")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$sugar_needles, col = "green")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$sugar_phloem, col = "brown")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$sugar_roots, col = "black")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$sugar_xylem_sh, col = "red")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$sugar_xylem_st, col = "orange")
    abline(h = 0.41, lty = 2, col = "blue")
    text(30, 0.43, "Expected Equilibrium", col = "blue", cex = 0.75)
    abline(h = SCb, lty = 2, col = "pink")
    text(25, SCb + 0.02, "\"bloom\" threshold", col = "pink", cex = 0.75)

    plot(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$starch_needles +
           sugar_model_cpp$Sugar$starch_phloem +
           sugar_model_cpp$Sugar$starch_roots +
           sugar_model_cpp$Sugar$starch_xylem_sh +
           sugar_model_cpp$Sugar$starch_xylem_st, main = "Starch", col = "blue", type = "l", ylim = c(0, 1.7), xlab = "Days of the Year", ylab = "Sugar, kg C")
    abline(h = 0, lty = 2, col = "grey")
    abline(v = 0, lty = 2, col = "grey")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$starch_phloem +
            sugar_model_cpp$Sugar$starch_roots +
            sugar_model_cpp$Sugar$starch_xylem_sh +
            sugar_model_cpp$Sugar$starch_xylem_st, col = "blue")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$starch_needles, col = "green")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$starch_phloem, col = "brown")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$starch_roots, col = "black")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$starch_xylem_sh, col = "red")
    lines(sugar_model_cpp$Growth$year, sugar_model_cpp$Sugar$starch_xylem_st, col = "orange")

  } else {
    CASSIA_sugar <- CASSIA_out[[1]]
    CASSIA_sugar$date <- as.POSIXct(as.character(CASSIA_sugar$date), format = "%Y-%m-%d")
    rownames(CASSIA_sugar) <- CASSIA_sugar$date

    ## TODO: add the data here!
    par(mfrow = c(2, 1))
    plot(CASSIA_sugar$date, CASSIA_sugar$sugar.needles + CASSIA_sugar$sugar.phloem + CASSIA_sugar$sugar.roots, main = "Sugar", col = "blue", type = "l", ylim = c(0, 0.8), xlab = "Days of the Year", ylab = "Sugar, kg C")
    abline(h = 0, lty = 2, col = "grey")
    lines(CASSIA_sugar$date, CASSIA_sugar$sugar.needles + CASSIA_sugar$sugar.phloem + CASSIA_sugar$sugar.roots, col = "blue")
    lines(CASSIA_sugar$date, CASSIA_sugar$sugar.needles, col = "green")
    lines(CASSIA_sugar$date, CASSIA_sugar$sugar.phloem, col = "brown")
    lines(CASSIA_sugar$date, CASSIA_sugar$sugar.roots, col = "black")
    abline(h = 0.41, lty = 2, col = "blue")
    text(30, 0.43, "Expected Equilibrium", col = "blue", cex = 0.75)
    abline(h = SCb, lty = 2, col = "pink")
    text(25, SCb + 0.02, "\"bloom\" threshold", col = "pink", cex = 0.75)

    plot(CASSIA_sugar$date, CASSIA_sugar$starch.needles + CASSIA_sugar$starch.phloem + CASSIA_sugar$starch.roots, main = "Starch", col = "blue", type = "l", ylim = c(0, 1.7), xlab = "Days of the Year", ylab = "Sugar, kg C")
    abline(h = 0, lty = 2, col = "grey")
    abline(v = 0, lty = 2, col = "grey")
    lines(CASSIA_sugar$date, CASSIA_sugar$starch.needles + CASSIA_sugar$starch.phloem + CASSIA_sugar$starch.roots, col = "blue")
    lines(CASSIA_sugar$date, CASSIA_sugar$starch.needles, col = "green")
    lines(CASSIA_sugar$date, CASSIA_sugar$starch.phloem, col = "brown")
    lines(CASSIA_sugar$date, CASSIA_sugar$starch.roots, col = "black")
  }
}
