library(testthat)
library(CASSIA)

# Shared fixture: 2015-2017 weather (1096 days, 2016 is a leap year)
data("weather_original")
weather_original$dates <- as.Date(
  strptime(paste(rep(2015:2017, times = c(365, 366, 365)), weather_original$X),
           format = "%Y %j"))

# ---------------------------------------------------------------------------
# Output structure — catches API breakage after refactoring
# ---------------------------------------------------------------------------

test_that("CASSIA_cpp returns expected list structure", {
  out <- CASSIA_cpp(weather = weather_original, site = "Hyde", tests = TRUE)
  expect_true(all(c("Growth", "Sugar", "Preles") %in% names(out)))
  expect_s3_class(out$Growth, "data.frame")
  expect_true(all(c("bud_growth", "diameter_growth", "needle_growth",
                    "root_growth", "height_growth",
                    "respiration_growth", "respiration_maintenance") %in%
                    names(out$Growth)))
  expect_equal(nrow(out$Growth), nrow(weather_original))
})

test_that("organ_level_sugar model Sugar output has all organ columns", {
  out <- CASSIA_cpp(weather = weather_original, site = "Hyde",
                    organ_level_sugar = TRUE, tests = TRUE)
  organs <- c("needles", "phloem", "roots", "xylem_sh", "xylem_st")
  expected <- c(paste0("sugar_",  organs), paste0("starch_", organs))
  expect_true(all(expected %in% names(out$Sugar)))
})

# ---------------------------------------------------------------------------
# Model logic — parameter-independent response direction tests.
# These test whether the model equations behave sensibly, not whether
# the parameters are calibrated correctly (that is what tests_c++_against_R
# is for, with its visual plots against Hyde_daily_original).
# ---------------------------------------------------------------------------

test_that("zeroing GPP substantially reduces total organ growth", {
  # The model draws on stored reserves even when current GPP = 0, so growth
  # does not collapse to exactly zero. The test checks that growth is at least
  # 90 % lower than the baseline run, not that it is precisely zero.
  weather_zero <- weather_original
  weather_zero$P <- 0

  total_growth <- function(out) {
    g <- out$Growth
    sum(g$needle_growth + g$diameter_growth + g$root_growth + g$height_growth)
  }

  run_base <- CASSIA_cpp(weather = weather_original, site = "Hyde", tests = TRUE)
  run_zero <- CASSIA_cpp(weather = weather_zero,    site = "Hyde", tests = TRUE)
  expect_lt(total_growth(run_zero), total_growth(run_base) * 0.1)
})

test_that("doubling GPP increases total organ growth", {
  weather_high <- weather_original
  weather_high$P <- weather_original$P * 2

  total_growth <- function(out) {
    g <- out$Growth
    sum(g$needle_growth + g$diameter_growth + g$root_growth + g$height_growth)
  }

  run_base <- CASSIA_cpp(weather = weather_original, site = "Hyde", tests = TRUE)
  run_high <- CASSIA_cpp(weather = weather_high,    site = "Hyde", tests = TRUE)
  expect_gt(total_growth(run_high), total_growth(run_base))
})

test_that("organ_level_sugar model: higher GPP produces higher mean needle sugar", {
  weather_high <- weather_original
  weather_high$P <- weather_original$P * 2

  run_base <- CASSIA_cpp(weather = weather_original, site = "Hyde",
                          organ_level_sugar = TRUE, tests = TRUE)
  run_high <- CASSIA_cpp(weather = weather_high, site = "Hyde",
                          organ_level_sugar = TRUE, tests = TRUE)

  expect_gt(mean(run_high$Sugar$sugar_needles, na.rm = TRUE),
            mean(run_base$Sugar$sugar_needles, na.rm = TRUE))
})
