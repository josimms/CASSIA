Cpp_code_test <- function() {
  ## TODO: again make sure that this makes sense in the context and also make sure that this is if hte sperling model is true!

  if (phloem.trigger == T) {
    if (count%%2 != 0) {
      if (length(which(sugar.phloem < sperling[c("SCb"),c(site)])) == 0) {
        warning(paste("Never cold enough for the model to trigger bud burst, sugar never lower than", sperling[c("SCb"),c(site)], "kg C"))
      } else {
        sB0 <- which(sugar.phloem < sperling[c("SCb"),c(site)])[1]
      }
    }
  } else {
    if (count%%2 != 0) {
      if (length(which(sugar.needles+sugar.phloem+sugar.roots+sugar.xylem.sh+sugar.xylem.st < sperling[c("SCb"),c(site)])) == 0) {
        warning(paste("Never cold enough for the model to trigger bud burst, sugar never lower than", sperling[c("SCb"),c(site)], "kg C"))
      } else {
        sB0 <- which(sugar.needles+sugar.phloem+sugar.roots+sugar.xylem.sh+sugar.xylem.st < sperling[c("SCb"),c(site)])[1]
      }
    }
  }

  ## TODO: make sure that the variable names are changed here, or at least make sense with the output from the sugar model
  if (storage.reset == FALSE) {
    if (n.year != 1) {
      sperling[c("starch.needles0"), c(site)] = starch.needles[n.days]
      sperling[c("starch.phloem0"), c(site)] = starch.phloem[n.days]
      sperling[c("starch.xylem.sh0"), c(site)] = starch.xylem.sh[n.days]
      sperling[c("starch.xylem.st0"), c(site)] = starch.xylem.st[n.days]
      sperling[c("starch.roots0"), c(site)] = starch.roots[n.days]

      sperling[c("sugar.needles0"), c(site)] = sugar.needles[n.days]
      sperling[c("sugar.phloem0"), c(site)] = sugar.phloem[n.days]
      sperling[c("sugar.xylem.sh0"), c(site)] = sugar.xylem.sh[n.days]
      sperling[c("sugar.xylem.st0"), c(site)] = sugar.xylem.st[n.days]
      sperling[c("sugar.roots0"), c(site)] = sugar.roots[n.days]

      sperling[c("Ad0.needles"),c(site)] = Ad.needles[n.days]
      sperling[c("Ad0.phloem"),c(site)] = Ad.phloem[n.days]
      sperling[c("Ad0.roots"),c(site)] = Ad.roots[n.days]
      sperling[c("Ad0.xylem.sh"),c(site)] = Ad.xylem.sh[n.days]
      sperling[c("Ad0.xylem.st"),c(site)] = Ad.xylem.st[n.days]
    }
  }
}


