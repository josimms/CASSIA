## CASSIA Model

### Ideas and structure

CASSIA model is an intra-annual growth model for an individual tree in boreal conditions. Seasonal organ level cell growth is modelled, as well as sugar and water when the appropriate settings are chosen. Further deatails for the indervidual functions can be found with the functions themselves.

The main structure and equations are found in Schiestl‐Aalto 2015 where the science behind this model as well as the basic principle and structure are clearly explained. The variable links in the papers and the model are written in the vingette section of this package.

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

### Example

As the package includes preprocessed weather data from Hyytiälä (via SPP model further information in vingettes) and amongst others the Hyytiälä configuration it is possible to simply run the model by stipulating these two arguments. This will run the model with its most basic functions, although additional functions can be easily added by toggles as seen in the second example. 

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
It is likely that the values you have chosen for the parameters have caused one of the outputs to not make sense. Thus the bounds of the parameters should be considered very carefully. If the problem persists, then report the error.

## Literature
### (please add if any missing!)

Ding, Yiyang, et al. "Temperature and moisture dependence of daily growth of Scots pine (Pinus sylvestris L.) roots in Southern Finland." Tree Physiology 40.2 (2020): 272-283.

Schiestl-Aalto, Pauliina, et al. "Analysis of the NSC storage dynamics in tree organs reveals the allocation to belowground symbionts in the framework of whole tree carbon balance." Frontiers in Forests and Global Change 2 (2019): 17.

Schiestl-Aalto, Pauliina, and Annikki Mäkelä. "Temperature dependence of needle and shoot elongation before bud break in Scots pine." Tree Physiology 37.3 (2017): 316-325.

Schiestl‐Aalto, Pauliina, et al. "CASSIA–a dynamic model for predicting intra‐annual sink demand and interannual growth variation in S cots pine." New Phytologist 206.2 (2015): 647-659.

Schiestl-Aalto, Pauliina, et al. "Physiological growth model CASSIA predicts carbon allocation and wood formation of Scots pine." CyberPlantS: a European initiative towards collaborative plant modeling (2013): 159.

## Used in
### (please add if any missing!)

Ding, Yiyang, et al. "Distinct patterns of below-and aboveground growth phenology and litter carbon inputs along a boreal site type gradient." Forest Ecology and Management 489 (2021): 119081.

Hellén, Heidi, et al. "Sesquiterpenes and oxygenated sesquiterpenes dominate the VOC (C 5–C 20) emissions of downy birches." Atmospheric Chemistry and Physics 21.10 (2021): 8045-8066.

Taipale, Ditte, et al. "The importance of accounting for enhanced emissions of monoterpenes from new Scots pine foliage in models-A Finnish case study." Atmospheric Environment: X 8 (2020): 100097.

Taipale, Ditte, et al. "Emissions of monoterpenes from new Scots pine foliage: dependency on season, stand age and location and importance for models." Biogeosciences Discussions (2020): 1-42.

Tian, Xianglin, et al. "Disaggregating the effects of nitrogen addition on gross primary production in a boreal Scots pine forest." Agricultural and Forest Meteorology 301 (2021): 108337.

## Related Literature
### (please add if any missing!)

Sperling, Or, et al. "Predicting bloom dates by temperature mediated kinetics of carbohydrate metabolism in deciduous trees." Agricultural and Forest Meteorology 276 (2019): 107643.
