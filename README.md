# CASSIA — Intra-Annual Tree Growth Model

CASSIA simulates daily carbon allocation and organ-level growth for an individual
Scots pine in boreal conditions. Seasonal growth of needles, shoots, diameter, and
roots is modelled, with optional sugar dynamics, xylogenesis, and mycorrhizal
interactions. Built-in parameters and weather data for Hyytiälä (SMEAR II) allow
the model to run out of the box.

> **Branches:** `main` is the stable, published version.
> `adding_externals` is the development branch with experimental features
> (enzyme-driven sugar model, p-hydro photosynthesis, mycorrhizal coupling).
> For support contact Joanna Simms: joanna.x.simms@helsinki.fi

---

## For Users

### Install

```r
install.packages("devtools")
devtools::install_github("josimms/CASSIA")
```

### Quick start

The package includes built-in Hyytiälä weather data (2015–2017), so you can run
the model immediately after installing:

```r
library(CASSIA)

# Add a date column to the built-in weather data
weather_original$dates <- as.Date(
  strptime(paste(rep(2015:2017, times = c(365, 366, 365)), weather_original$X),
           format = "%Y %j"))

# Run the model with default settings
out <- CASSIA_cpp(weather = weather_original, site = "Hyde")
```

The output is a list with three dataframes:

| Element | Contents |
|---------|----------|
| `out$Growth` | Daily carbon fluxes and organ growth (bud, wall, needle, root, height; respiration; storage; mycorrhiza) |
| `out$Sugar` | Daily sugar and starch pool values per organ (populated when the sugar model is active) |
| `out$Preles` | Daily photosynthesis outputs: GPP, ET, soil water |

```r
head(out$Growth)
```

See `?CASSIA_cpp` for the full list of arguments and what they do.
See `vignette("Running_The_Model")` for a guided walkthrough.

### Changing settings (toggles)

Model behaviour is controlled by logical toggle arguments. For example, to turn
off the GPP correction on needle elongation:

```r
out2 <- CASSIA_cpp(weather = weather_original, site = "Hyde", LN_estim = FALSE)
```

The full list of toggles is in `?CASSIA_cpp`. If you choose a combination of
toggles that cannot coexist, the model will warn you and adjust automatically —
check those warnings carefully.

### Changing parameters

Default parameters are stored in built-in dataframes (`parameters_p`, `ratios_p`,
`common_p`, `sperling_p`, `repo_p`). Copy and modify the one you need:

```r
parameters_updated <- parameters_p
parameters_updated["root.lifetime", "Hyde"] <- 2
out3 <- CASSIA_cpp(weather = weather_original, site = "Hyde",
                   parameters = parameters_updated)
```

All parameters are explained in the CASSIA instruction booklet
(`CASSIA_Instruction_Booklet.pdf`, in the package root or on GitHub) and in
the original articles listed below.

### Using your own weather data

The weather dataframe must have these columns (in these units):

| Column | Description | Units |
|--------|-------------|-------|
| `X` | Day of year | integer |
| `T` | Air temperature | °C |
| `P` | Photosynthesis per tree (when `photosynthesis_as_input = TRUE`) | g C m⁻² day⁻¹ |
| `TSA` | Soil temperature, A horizon | °C |
| `TSB` | Soil temperature, B horizon | °C |
| `MB` | Soil moisture, B horizon | m³ m⁻³ |
| `Rain` | Precipitation | mm day⁻¹ |
| `dates` | Calendar date | `Date` class |

The scripts `R/Hyytiala_Weather_Fast.R` and `R/ERAS_weather_processing.R` show
how raw station and ERA5 data are formatted. The model will give an informative
error if column names are wrong or NAs are present.

### Calibrating for a new site

If you are not using Hyytiälä (`"Hyde"`) or Lettosuo, the model is not yet
calibrated for your site. See `vignette("Sensitivity_Analysis")` for guidance on
sensitivity analysis and parameter fitting.

### Vignettes

| Vignette | Topic |
|----------|-------|
| `vignette("Running_The_Model")` | Full walkthrough: install, run, interpret output, change parameters |
| `vignette("Visualisation")` | Plotting model output and comparing to observations |
| `vignette("Datasets")` | Exploring and modifying built-in parameter tables and datasets |
| `vignette("Component_Models")` | Using Repola allometry and PRELES photosynthesis independently |
| `vignette("Weather")` | Preparing and formatting weather data for CASSIA |
| `vignette("ERA5_Weather")` | ERA5 reanalysis workflow and the `calculate_VPD` function |
| `vignette("ring_width_process")` | Step-by-step walkthrough of the xylogenesis sub-model |
| `vignette("Sensitivity_Analysis")` | Testing parameter sensitivity and calibrating for a new site |

---

## For Developers and Contributors

### Package layout

```
CASSIA/
├── R/                    # Exported R functions + roxygen2 documentation blocks
├── src/                  # C++ model core, compiled via Rcpp
├── inst/include/         # C++ header files shared across translation units
├── data/                 # Built-in datasets (weather_original, parameters_p, …)
├── data-raw/             # Scripts that built the data/ objects
├── man/                  # Auto-generated help files — do not edit by hand
├── vignettes/            # User-facing tutorials (rendered on CRAN/GitHub)
├── external/             # External dependency code
└── development_notebooks/ # Developer notebooks (model comparisons, calibration)
```

The C++ model (`src/`) is compiled by Rcpp on install. The R wrapper
`R/CASSIA_cpp.R` handles argument validation and calls the compiled function.

### Regenerating documentation

Documentation for the main user-facing functions is written as roxygen2 `#'`
blocks in the `R/` files and auto-generated into `man/`. After editing any `#'`
block, regenerate with:

```r
devtools::document()
```

### Adding a new site

1. Add a column to `parameters_p`, `ratios_p`, and `common_p` with the site name.
2. Rebuild the data objects:
   ```r
   usethis::use_data(parameters_p, overwrite = TRUE)
   usethis::use_data(ratios_p, overwrite = TRUE)
   usethis::use_data(common_p, overwrite = TRUE)
   ```
3. Pass the new site name to `CASSIA_cpp(..., site = "YourSite")`.

The `data-raw/` scripts show how the built-in parameter tables were constructed —
use them as templates.

### Adding a new toggle

1. Add the argument to `CASSIA_cpp()` in `R/CASSIA_cpp.R` with a default of `FALSE`.
2. Add a `#' @param your_toggle Logical. ...` roxygen2 line to the documentation block.
3. Pass the argument through to the C++ function in `src/` and implement the behaviour.
4. Run `devtools::document()` to update `man/CASSIA_cpp.Rd`.

### C++/R interface

The model core is a single Rcpp function defined in `src/`. It returns a named
list:

```cpp
return Rcpp::List::create(Rcpp::_["Growth"] = df,
                           Rcpp::_["Sugar"]  = df2,
                           Rcpp::_["Preles"] = df3);
```

Access these on the R side as `out$Growth`, `out$Sugar`, `out$Preles`.

### Ongoing projects

- **Lettosuo parameterisation** — Alexis Lehtonen
- **Enzyme-driven sugar model + mycorrhizal coupling** — Joanna Simms

---

## Literature

Schiestl-Aalto, P., et al. (2015). CASSIA — a dynamic model for predicting
intra-annual sink demand and interannual growth variation in Scots pine.
*New Phytologist* 206(2): 647–659.

Schiestl-Aalto, P., and Mäkelä, A. (2017). Temperature dependence of needle and
shoot elongation before bud break in Scots pine. *Tree Physiology* 37(3): 316–325.

Schiestl-Aalto, P., et al. (2019). Analysis of the NSC storage dynamics in tree
organs reveals the allocation to belowground symbionts in the framework of whole
tree carbon balance. *Frontiers in Forests and Global Change* 2: 17.

Ding, Y., et al. (2020). Temperature and moisture dependence of daily growth of
Scots pine roots in Southern Finland. *Tree Physiology* 40(2): 272–283.

## Used in

Ding, Y., et al. (2021). Distinct patterns of below- and aboveground growth
phenology and litter carbon inputs along a boreal site type gradient.
*Forest Ecology and Management* 489: 119081.

Tian, X., et al. (2021). Disaggregating the effects of nitrogen addition on gross
primary production in a boreal Scots pine forest. *Agricultural and Forest
Meteorology* 301: 108337.

Taipale, D., et al. (2020). The importance of accounting for enhanced emissions of
monoterpenes from new Scots pine foliage in models — a Finnish case study.
*Atmospheric Environment: X* 8: 100097.
