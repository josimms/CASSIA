## CASSIA Model

### Ideas and structure

CASSIA model is an intra-annual growth model for an individual tree in boreal conditions. Seasonal organ level cell growth is modelled, as well as sugar and water when the appropriate settings are chosen. Further deatails for the indervidual functions can be found with the functions themselves.

The main mathematical structure and equations are found in Schiestl‐Aalto 2015 where the science behind this model as well as the basic principle and structure are clearly explained. The variable links in the papers and the model are written in the vingette section of this package. This basic equations have been added to and reported in later publications listed below. This package also has newer developments not yet published in papers such as a sugar internal allocation model and xylogenesis.

This model has been used in numerous papers (see Literature) - mainly considering Hyytiälä (SMEAR II Station, University of Helsinki). Most of the code has been written by Schiestl‐Aalto, with only small additions by others. The development of the code can be seen in the papers below, although this package does have newer developments not yet published in papers such as a sugar addition and xylogenesis. Note: future developments of the CASSIA model will seperate these into spereate functions, but at the moment all of the code is together in one file to be able to make the original code into a package.

### Ongoing Projects

# Alexis Lehtonen
Parameterisation for the CASSIA Lehtosuo site

# Joanna Simms 
Addition of a enzyme driven sugar model and mycorrhizal interactions. Thus linking with a soil model (SYMPHONY) and mycorrhizal model (MYCOFON).

![kuva](https://github.com/josimms/MycoModel/assets/102613042/1a465070-6995-4f73-bef7-4e7920bca289)

The full details should be added when the code is finished, but indervidual functions are written in the function help files.

### Downloading Package for Use

```{r}
install.packages("devtools")  # source: https://github.com/hoxo-m/githubinstall
library(devtools)
install_github("josimms/CASSIA")
```

### Example R

As the package includes preprocessed weather data from Hyytiälä (via SPP model further information in vingettes) and amongst others the Hyytiälä configuration it is possible to simply run the model by stipulating these two arguments. This will run the model with its most basic functions, although additional functions can be easily added by toggles as seen in the second example. Toggles are found in the documentation.

```{r}
library(CASSIA)
CASSIA(weather = Hyde_weather, site = "Hyde")
```

Hyde_weather is included in the package, and another function is included Hyde_Data_Creation, which could be used as a reference for building your own weather data.

Here the argument for mychorrhiza has been changed to FALSE (default is TRUE). Now mychorrhiza settings will be taken into account when the model is running. To see all of the possible processes that are considered in the model look to the CASSIA function documentation. This also provides the initial conditions and settings of the model. The formatting and inputs are also listed here and should be considered.

```{r}
library(CASSIA)
CASSIA(weather = Hyde_weather, site = "Hyde", mychorrhiza = FALSE)
```

If you have an error along the lines of 

```{r}
Error in if GENERIC ARGUMENT missing value where TRUE/FALSE needed
```
It is likely that the values you have chosen for the parameters have caused one of the outputs to not make sense. Thus the bounds of the parameters should be considered very carefully. If the problem persists, then report then send joanna.x.simms@helsinki.fi an email.

(Code working as of 19th Jan 2024 - contact Joanna if not working)

### Example C++ model via R interface

The C++ model has less automatic features than the R version of the model. This means that when you call the function you have to be more explicit about all of the arguments as well as including different weather data. There is a working example in the package, however this is not fully documented. A basic example is as follows.

(Code currently under revision, so could not be working - this code will be updated by Feburary with a working version)

```{r}
### Toggle setting
storage_rest = T
storage_grows = F
LH_estim = T
LN_estim = T
mN_varies = T
LD_estim = T
sD_estim_T_count = F
trees_grow = F
growth_decreases = F
needle_mass_grows = F
mycorrhiza = T
root_as_Ding = T
sperling_sugar_model = F
using_spp_photosynthesis = T
xylogensis_option = F
environmental_effect_xylogenesis = F
temp_rise = F
drought = F
Rm_acclimation = F
etmodel = F
LOGFLAG = F

### Non automatic parameters
N_parameters = c(1, 1)
pPREL = c(413.0, 0.450, 0.118, 3.0, 0.748464, 12.74915, -3.566967, 18.4513, -0.136732,
            0.033942, 0.448975, 0.500, -0.364, 0.33271, 0.857291, 0.041781,
            0.474173, 0.278332, 1.5, 0.33, 4.824704, 0.0, 0.0, 180.0, 0.0, 0.0, 10.0,
            -999.9, -999.9, -999.9)

### Weather dataset updated with the extra terms needed for the C++ model
weather_original_2015 = read.csv(file = "./data/weather_original_2015.csv", header = T, sep = ",")
weather_original_2016 = read.csv(file = "./data/weather_original_2016.csv", header = T, sep = ",")
weather_original_2017 = read.csv(file = "./data/weather_original_2017.csv", header = T, sep = ",")
weather_original = rbind(rbind(weather_original_2015, weather_original_2016), weather_original_2017)

extras = data.frame(Nitrogen = rep(0.012, length = nrow(weather_original)),
                    PAR = data_format[substring(data_format$Date, 1, 4) %in% 2015:2017,c("PAR")],
                    VPD = data_format[substring(data_format$Date, 1, 4) %in% 2015:2017,c("VPD")],
                    CO2 = data_format[substring(data_format$Date, 1, 4) %in% 2015:2017,c("CO2")],
                    fAPAR = rep(0.7, length = nrow(weather_original)))
weather_original <- cbind(weather_original, extras)
weather_original <- weather_original[-c(365+365),]

### Function call
CASSIA_yearly(2015, 2016, weather_original, GPP_ref,
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
```

The other parameters are defined in the model automatically, although they have to be called, and can be found in the data folder.

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

