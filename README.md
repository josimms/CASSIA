## CASSIA Model

NOTE: The CASSIA+ Instruction and Plans Booklet has the latest information in terms of what the model will include and what the model currently includes. If you have any questions please email: joanna.x.simms@helsinki.fi

### Ideas and structure

CASSIA model is an intra-annual growth model for an individual tree in boreal conditions. Seasonal organ level cell growth is modelled, as well as sugar and water when the appropriate settings are chosen. Further deatails for the indervidual functions can be found with the functions themselves.

The main mathematical structure and equations are found in Schiestl‐Aalto 2015 where the science behind this model as well as the basic principle and structure are clearly explained. The variable links in the papers and the model are written in the vingette section of this package. This basic equations have been added to and reported in later publications listed below. This package also has newer developments not yet published in papers such as a sugar internal allocation model and xylogenesis.

This model has been used in numerous papers (see Literature) - mainly considering Hyytiälä (SMEAR II Station, University of Helsinki). Most of the code has been written by Schiestl‐Aalto, with small additions by others. Note: future developments of the CASSIA model will seperate different versions of the model into different functions.

### Downloading Package for Use

```{r}
install.packages("devtools")  # source: https://github.com/hoxo-m/githubinstall
library(devtools)
install_github("josimms/CASSIA")
```

### Example R

**markdown files cover using the model, also sensitivity analysis and weather processing**

The CASSIA model (function: CASSIA_cpp) can be simply run with two arguments, as the package includes preprocessed parameters and weather data from Hyytiälä (via SPP model further information in vingettes). This will run the model with its most basic functions, although additional functions can be easily added by toggles as seen in the second example. Toggles are found in the documentation. 

Some weather data is included in the package and formated by "process_weather_data". Other functions such as Hyde_Data_Creation and raw_to_daily_monthly_hyytiala, are used to get the weather data in the correct form. These can also be used as a reference for building your own weather data. NOTE: currently the weather processes are being rewritten for the new photosynthesis functionalities. 

Wather data format also depends on the photsynthesis model used, the basic input of the model is that photosynthesis per tree is an input, but it can also be calculated in the CASSIA package if needed. The different photosynthesis input methods need different weather data formats. The weather data needed is detailed in the represctive photosynthesis models.
1. PRELES
2. p-hydro
The R function will give an error if the weather data is wrongly named, within unexpected bounds or if it has NAs.

```{r}
# import the package
library(CASSIA)

weather_original <- read_and_combine_weather_data(2015, 2018, base_path = "~/Documents/CASSIA/data/weather_original_")
# TODO: change the path to your path. Data is found in ./data in the package
weather_original$date <- as.Date(strptime(paste(rep(2015:2018, times = c(365, 366, 365, 365)), weather_original$X), format = "%Y %j"))

# Run the model
out = CASSIA_cpp(weather = weather_original, site = "Hyde")
```

To use a different model setting you can use a toggle. An example that doesn't require changing any of the rest of the inputs is LN.estim. This means that a GPP correction is not applied to the potential needle growth calculations.

```{r}
CASSIA_cpp(weather = weather_original, site = "Hyde", LN.estim = FALSE)
```
There should be an error or a warning if you have chosen a set of toggles that cannot coexist. Usually the model will assume which combination you ment and run with the corrected toggles, so these should be checked carefully.

The package automatically has tables of parameters included. These are notated with a _p after the argument name as in the function. The easiest way to reformulate the model is thus take these as a base and change the parameters you would like to accordingly. All of the parameters are named in the code objects and then are explained further in the CASSIA instruction booklet and the original articles.

In the data-raw folder you can also see how the _p objects are made if you would prefer this as a basis for your new parameters. An example of the first method is:

```{r}
ratios_new = ratios_p
ratios_new[1,c("Hyde")] = 0.65
CASSIA_cpp(weather = weather_original, site = "Hyde", ratios = ratios_new)
```

If you are using the model in Hyytiälä or Väriö then the model is calibrated. If not the model is not yet calibrated for your site. This means that you need to calibrate the model! To do this look at the markdown files for advice.

**markdown files cover using the model again and sensitivity analysis**

Sidenote: If you have an error along the lines of 

```{r}
Error in if GENERIC ARGUMENT missing value where TRUE/FALSE needed
```
It is likely that the values you have chosen for the parameters have caused one of the outputs to not make sense. Thus the bounds of the parameters should be considered very carefully. If the problem persists, then report then send joanna.x.simms@helsinki.fi an email.

(Code working as of 19th Nov 2024 - contact Joanna if not working)

### Weather Data

It is hard to make an example file for the weather data needed as all weather data sources are different. Examples of such processes can be found in the files Hyytiala_Weather_Fast.R and ERAS_weather_processing.R - although the data needs to be downloaded by the user first.

The basic process for all of the data is to:
1. Download the data.
2. Format the data as in the function help files in the units required (?CASSIA).
3. Gapfill / Bias correct. If the data is taken from local weather stations you will need to gapfill the results with the closely stations. If you have taken the result from a global data set you will then need to bias correct with the local sites. (Gapfill example in Hyyiälä data porcessing and the bias correction in the ERA5 file)

More information can be found in Running_The_Model.Rmd

### Ongoing Projects

#### Parameterisation for the CASSIA Lehtosuo site
Alexis Lehtonen

#### Addition of a enzyme driven sugar model and mycorrhizal interactions. Thus linking with a soil model (SYMPHONY) and mycorrhizal model (MYCOFON). 
Joanna Simms

![kuva](https://github.com/josimms/MycoModel/assets/102613042/1a465070-6995-4f73-bef7-4e7920bca289)

The full details should be added when the code is finished, but indervidual functions are included in the package with relevent help files.


## Literature
Ding, Yiyang, et al. "Temperature and moisture dependence of daily growth of Scots pine (Pinus sylvestris L.) roots in Southern Finland." Tree Physiology 40.2 (2020): 272-283.

Schiestl-Aalto, Pauliina, et al. "Analysis of the NSC storage dynamics in tree organs reveals the allocation to belowground symbionts in the framework of whole tree carbon balance." Frontiers in Forests and Global Change 2 (2019): 17.

Schiestl-Aalto, Pauliina, and Annikki Mäkelä. "Temperature dependence of needle and shoot elongation before bud break in Scots pine." Tree Physiology 37.3 (2017): 316-325.

Schiestl‐Aalto, Pauliina, et al. "CASSIA–a dynamic model for predicting intra‐annual sink demand and interannual growth variation in S cots pine." New Phytologist 206.2 (2015): 647-659.

Schiestl-Aalto, Pauliina, et al. "Physiological growth model CASSIA predicts carbon allocation and wood formation of Scots pine." CyberPlantS: a European initiative towards collaborative plant modeling (2013): 159.

## Used in
Ding, Yiyang, et al. "Distinct patterns of below-and aboveground growth phenology and litter carbon inputs along a boreal site type gradient." Forest Ecology and Management 489 (2021): 119081.

Hellén, Heidi, et al. "Sesquiterpenes and oxygenated sesquiterpenes dominate the VOC (C 5–C 20) emissions of downy birches." Atmospheric Chemistry and Physics 21.10 (2021): 8045-8066.

Taipale, Ditte, et al. "The importance of accounting for enhanced emissions of monoterpenes from new Scots pine foliage in models-A Finnish case study." Atmospheric Environment: X 8 (2020): 100097.

Taipale, Ditte, et al. "Emissions of monoterpenes from new Scots pine foliage: dependency on season, stand age and location and importance for models." Biogeosciences Discussions (2020): 1-42.

Tian, Xianglin, et al. "Disaggregating the effects of nitrogen addition on gross primary production in a boreal Scots pine forest." Agricultural and Forest Meteorology 301 (2021): 108337.

## Related Literature
Meyer, Astrid, et al. "Simulating mycorrhiza contribution to forest C-and N cycling-the MYCOFON model." Plant and soil 327 (2010): 493-517.

Perveen, Nazia, et al. "Priming effect and microbial diversity in ecosystem functioning and response to global change: a modeling approach using the SYMPHONY model." Global change biology 20.4 (2014): 1174-1190.

Sperling, Or, et al. "Predicting bloom dates by temperature mediated kinetics of carbohydrate metabolism in deciduous trees." Agricultural and Forest Meteorology 276 (2019): 107643.

