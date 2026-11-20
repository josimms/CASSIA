# Step-wise run (cassia_init / cassia_step / cassia_outputs): must give exactly the batch result,
# and runs must not share hidden state.

library(testthat)
library(CASSIA)

data("weather_original")
weather_original$dates <- as.Date(
  strptime(paste(rep(2015:2017, times = c(365, 366, 365)), weather_original$X),
           format = "%Y %j"))

settings <- list(
  default    = list(),
  preles     = list(photosynthesis_as_input = FALSE, preles = TRUE,
                    storage_reset = FALSE, soil_water_reset = FALSE),
  pipe_model = list(pipe_model_correction = TRUE, root_turnover = TRUE)
)

batch <- function(s, weather = weather_original)
  suppressWarnings(do.call(CASSIA_cpp, c(list(weather = weather, site = "Hyde", CASSIA_graphs = FALSE), s)))
start <- function(s, weather = weather_original)
  do.call(cassia_init, c(list(weather = weather, site = "Hyde", CASSIA_graphs = FALSE), s))
quiet <- function(expr) { capture.output(r <- suppressWarnings(expr)); r }

test_that("stepping in chunks to the end equals the batch run exactly", {
  for (name in names(settings)) {
    run <- start(settings[[name]])
    for (n in c(1, 17, 365)) quiet(cassia_step(run, n))
    quiet(cassia_step(run, nrow(weather_original)))
    expect_true(cassia_day(run)$done)
    expect_identical(cassia_outputs(run), quiet(batch(settings[[name]])), info = name)
  }
})

test_that("outputs part-way are the first days of the full run", {
  run <- start(settings$preles)
  quiet(cassia_step(run, 400))
  part <- cassia_outputs(run)
  full <- quiet(batch(settings$preles))
  expect_equal(nrow(part$Growth), 400)
  expect_equal(cassia_day(run)$date, weather_original$dates[401])
  for (df in names(full)) {
    expect_identical(part[[df]], {x <- full[[df]][1:400, , drop = FALSE]; rownames(x) <- NULL; x}, info = df)
  }
})

test_that("two runs stepped side by side equal the same runs done separately", {
  warmer <- weather_original
  warmer$T <- warmer$T + 2
  a <- start(settings$preles)
  b <- start(settings$preles, warmer)
  while (!cassia_day(a)$done) {
    quiet(cassia_step(a, 7))
    quiet(cassia_step(b, 7))
  }
  quiet(cassia_step(b, nrow(warmer)))
  expect_identical(cassia_outputs(a), quiet(batch(settings$preles)))
  expect_identical(cassia_outputs(b), quiet(batch(settings$preles, warmer)))
})

test_that("cassia_set_forcing replaces the next day's drivers", {
  run <- start(settings$default)
  quiet(cassia_step(run, 100))
  same <- weather_original[101, c("TSA", "MB")]
  cassia_set_forcing(run, list(TSA = same$TSA, MB = same$MB))       # unchanged values
  quiet(cassia_step(run, nrow(weather_original)))
  expect_identical(cassia_outputs(run), quiet(batch(settings$default)))

  changed <- weather_original
  changed$TSA[101] <- changed$TSA[101] + 5
  run <- start(settings$default)
  quiet(cassia_step(run, 100))
  cassia_set_forcing(run, c(TSA = changed$TSA[101]))
  quiet(cassia_step(run, nrow(weather_original)))
  expect_identical(cassia_outputs(run), quiet(batch(settings$default, changed)))

  expect_error(cassia_set_forcing(run, list(TSA = 1)), "finished")
  run <- start(settings$default)
  expect_error(cassia_set_forcing(run, list(nonsense = 1)), "unknown driver")
})
