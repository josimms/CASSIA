\name{CASSIA}
\alias{CASSIA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CASSIA growth model
%%  ~~function to do ... ~~
}
\description{
This function runs CASSIA model.
}
\usage{
CASSIA(weather, site, ratios = ratios_p, parameters = parameters_p, 
       common = common_p, sperling = sperling_p, sperling = sperling_p, repo = repo_p,
      storage.reset = TRUE, storage.grows = FALSE,
      LN.estim = TRUE, mN.varies = TRUE,
      LD.estim = TRUE, sd.estim.T.count = FALSE,
      LH.estim = TRUE,
      trees_grow = FALSE, growth_decreases = FALSE, needle_mass_grows = FALSE,
     phloem.trigger = FALSE,
      mychorrhiza = TRUE, root_as_Ding = TRUE,
      sperling_model = FALSE, xylogenesis = FALSE,
      PRELES_GOO = FALSE, environment_effect_xylogenesis = FALSE,
      photoparameters = 3, temp_rise =  FALSE, frought = FALSE, Rm_acclimation = TRUE,
      s.D0 = 79, s.H0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weather}{
    % TODO: was this mean or sums?
    
  This should be a daily dataframe input:
    Date: Year-Month-Day form
    Air Temoperature (T): Degree C
    Photosynthesis (P): g C m-2 day-1
    Temperature of the soil at soil horizon A (TSA): Degree C
    Temperature of the soil at soil horizon B (TSB): Degree C
    Soil moisture (MB): In B horizon m3 m-3
    Precipitation (Rain): mm day-1
    
%%     ~~Describe \code{nYears} here~~
}
  \item{site}{
  As CASSIA is very site and species spercific parameterisations are for spercific sites. 
    The possible build in parameterisations are "Hyde" for Hyytiälä and "Lettosuo" for Lettosuo. 
      More sites can be added by adding a column in the parameters_p dataframe. (see parameters_p vingette)
    A future plan is to sperate site type and species.
%%     ~~Describe \code{pCROBAS} here~~
}
  \item{ratios}{
  In this dataframe of input parameters the ratios of growth are outlined (see vingette)
%%     ~~Describe \code{pPRELES} here~~
}
  \item{parameters}{
  In this dataframe of input parameters the indervidual site parameters are outlined (see vingette)
%%     ~~Describe \code{pPRELES} here~~
}
  \item{common}{
  In this dataframe common growth parameters between sites are outlined (see vingette)
%%     ~~Describe \code{pYASSO} here~~
}
  \item{sperling}{
  In this dataframe the sperling sugar model parameters as well as the original sugar parameters are outlined (see vingette)
%%     ~~Describe \code{pYASSO} here~~
}
  \item{repo}{
  In this dataframe the repola growth model parameters are outlined (see vingette)
%%     ~~Describe \code{siteInfo} here~~
}
  \item{storage.reset}{
  If this is TRUE, the sugar and starch stores reset to their original values every year.
  If this is FALSE the inital values each year are taken from the last storage value from the last simulated year.
    The initial value in both cases is given in the sperling parameters.
%%     ~~Describe \code{initVar} here~~
}
  \item{storage.grows}{
  If this is TRUE, the potential storage space of sugar and starch is updated related to the height and diameter equations.
  If this is FALSE, the potential storage stays at the same level inspite of growth.
%%     ~~Describe \code{thinning} here~~
}
  \item{LN.estim}{
  If this is TRUE, causes the LN to be calucated from previous years GPP data:
    LN <- parameters[c("LN0"), c(site)] * GPP_previous_sum[which(GPP_previous_sum[,1] == year), 2] / mean(GPP_previous_sum[,2])
  If this is FALSE, causes the LN to be calucated from coefficients:
    LN <- LN * growth_photo_coef
    Currently FALSE if xylogenesis is TRUE
%%     ~~Describe \code{initClearcut} here~~
}
  \item{mN.varies}{
    % TODO: how to make the patches of code look good?
  If this is TRUE, then repola model used in the calulation of RmN.a. 
  If this is FALSE, then repola model not used in the calulation of RmN.a.
%%     ~~Describe \code{PAR} here~~
}
  \item{LD.estim}{
  LD is calulated with diamater and hieght growth coefficients (diameter_growth_coefficient, height_growth_coefficient see parameters vingette)
  If this is TRUE, then LD is altered by the GPP of the preceding five days.
  If this is FALSE, then LD is non altered.
%%     ~~Describe \code{TAir} here~~
}
  \item{sD.estim.T.count}{
    % TODO: J could understand this better
  If this is TRUE, sD is simply cumsum(sDA > 0), but sDA changes to sDA = sD0.T.count + cumsum(g.sD.T) if xylogenesis = TRUE
  If this is FALSE, sD is calucated by parameters[c("sD0.Trad"), c(site)] + cumsum(g.sD.T)
%%     ~~Describe \code{VPD} here~~
}
  \item{LH.estim}{
  LH is calulated with diamater and hieght growth coefficients (diameter_growth_coefficient, height_growth_coefficient see parameters vingette)
  If this is TRUE, then LH is altered by the GPP of the preceding five days.
  If this is FALSE, then LH is non altered.
%%     ~~Describe \code{Precip} here~~
}
  \item{trees_grow}{
  If this is TRUE, D0 and h0 are updated in the parameters dataframe (see vingettes)
  If this is FALSE, D0 and h0 not updated and are considered constant for the entire time of simulation.
    Currently TRUE if xylogenesis is TRUE
%%     ~~Describe \code{CO2} here~~
}
  \item{growth_decreases}{
  If this is TRUE, height_growth_coefficient and diameter_growth_coefficient (ratios) are changed depedning on a linear apprximation given between the minimum and maximum values and the years 1997 and 2020.
  If this is FALSE, height_growth_coefficient and diameter_growth_coefficient mean value is repreated for the years 1997:2020.
%%     ~~Describe \code{P0} here~~
}
  \item{needle_mass_grows}{
   If this is TRUE, the needle mass is constant. Although calucated with the repola equations on the first interation.
   If this is FALSE, the repola model is used to calucate the needle mass.
%%     ~~Describe \code{soilC} here~~
}
  \item{phloem.trigger}{
  If this is TRUE, the bud burst / start of bud burst is calculated to start when sugar phloem concentrations first decrease to a set sugar leel (sperling parameters). If this value is never reached then the SCb parameter (sperling) is used instead and a warning is given.
  If this is FALSE, the bud burst is considered to be the day the sperling parameter says.
    Currently FALSE if sperling_model is FALSE as sperling_model needed for the more dynamic sugar levels.
%%     ~~Describe \code{weatherYasso} here~~
}
  \item{mychorrhiza}{
  If this is TRUE, parameters[c("growth.myco"),c(site)] holds.
  If this is FALSE, parameters[c("growth.myco"),c(site)] = 0
    Currenty FALSE if the sperling_model is TRUE as mycorrhiza is included explicitly in this model.
%%     ~~Describe \code{litterSize} here~~
}
  \item{root_as_Ding}{
  If this is TRUE, consinders the growth of roots with the ratos of fibrous and non-fibrous roots as well as explicit stages of development.
  If this is FALSE, considers only the temperature and moisture driving root growth.
%%     ~~Describe \code{soilCtot} here~~
}
  \item{sperling_model}{
  If this is TRUE, sugars are considered to be controlled by a enzyme driven sugar model explained in Splering, 2019 (References). Organ concentrations are also considered seperately with a concentration driven movement of sugars betewen pools.
  If this is FALSE, sugars are considered in the way outlined in Schiestl-Aalto 2019 (References) with a set equilibrium point and no organ seperation.
%%     ~~Describe \code{defaultThin} here~~
}
  \item{xylogenesis}{
  If this is TRUE, xylogenesis controls the cell growth with clear impacts on the needle cohorts, diameter and roots. Considers early and late growth.
  If this is FALSE, growth is controlled by the simpler functions in Schiestl-Aalto 2015 (References).
%%     ~~Describe \code{ClCut} here~~
}
  \item{PRELES_GPP}{
  If this is TRUE, then PRELES_GPP is used to calculate the growth_photo_coef
  If this is FALSE, then growth_photo_coef = 1
%%     ~~Describe \code{inDclct} here~~
}
  \item{environment_effect_xylogenesis}{
    % TODO: explain this better
  If this is TRUE, then there is an effect of tempreature in the formation of cells.
  If this is FALSE, then development depends on phases and preset parameters only.
%%     ~~Describe \code{inAclct} here~~
}
  \item{photoparameters}{
  Value. Is an input for the PRELES_GPP function.
%%     ~~Describe \code{yassoRun} here~~
}
  \item{temp_rise}{
  If this is TRUE, and Rm_acclimation is TRUE then the Rm_accl is modified from 1 to 0.85.
  If this is FALSE, RM_accl = 1.
%%     ~~Describe \code{yassoRun} here~~
}
  \item{drought}{
  Currently has no function.
%%     ~~Describe \code{yassoRun} here~~
}
  \item{Rm_acclimation}{
  If this is TRUE, and temp_rise is TRUE then the Rm_accl is modified from 1 to 0.85.
  If this is FALSE, RM_accl = 1.
%%     ~~Describe \code{yassoRun} here~~
}
  \item{s.D0}{
  Value. The day of the year that is used to start the temperature sum. Valid for Finland.
%%     ~~Describe \code{yassoRun} here~~
}
  \item{s.H0}{
  Value. The day of the year that is used to start shoot growth.
%%     ~~Describe \code{yassoRun} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{output}{
  A list of two outputs. First list item is an array with the daily values. Second is an arrany with the yearly values (culmatively summed).

  Expored yearly variables: \cr
  1."year" \cr
  2."total starch" summed for the sperling model, orginal parameter when sperling model isn't used \cr
  3."total sugar" summed for the sperling model, orginal parameter when sperling model isn't used \cr
  4."wall.tot" kg C \cr
  5."height.tot" kg C \cr
  6."needle.tot" kg C \cr
  7."tot.Rm"  kg C \cr
  8."tot.Rg" kg C \cr
  9."tot.P" sum g C m-2 day-1 Photosynthesis for one tree\cr
  10."cumsum.PF" sum of g C m-2 day-1 Photosynthesis for whole stand \cr
  11."cum.Daily.H.tot" mm \cr
  12."tot.mm" mm \cr
  13."needle_mass" C kg \cr
  14."sum.needle.cohorts" \cr
  15."sugar.needles" kg C\cr
  16."sugar.phloem" kg C\cr
  17."sugar.xylem.sh" kg C\cr
  18"sugar.xylem.st" kg C\cr
  19."sugar.roots" kg C\cr
  20."starch.needles" kg C\cr
  21."starch.phloem" kg C\cr
  22."starch.xylem.sh" kg C\cr
  23."starch.xylem.st" kg C\cr
  24."starch.roots" kg C)
    
  Expored daily variables: \cr
  1."date" Year-Month-Day \cr
  2."year" \cr
  3."day" \cr
  4."bud.tot.growth"  kg C \cr
  5."wall.tot.growth"  kg C \cr
  6."needle.tot.growth"  kg C \cr
  7."root.tot.growth" kg C \cr
  8."height.tot.growth" kg C \cr
  9."Rg.tot" kg C \cr
  10."Rm.tot" kg C \cr
  11."height.tot"  kg C \cr
  12."wall.tot"  kg C \cr
  13."storage"  kg C \cr
  14."sugar"  kg C \cr
  15."starch"  kg C \cr
  16."storage_term" % based on the maximum storage equation \cr
  17."to.mycorrhiza" kg C \cr
  18."mycorrhiza.tot" kg C \cr
  19."P" g C m-2 day-1 Photosynthesis for one tree \cr
  20."to_sugar"  kg C \cr
  21."to_starch"  kg C \cr
  22."Daily.H.tot" mm \cr
  23."Daily.N.tot" mm \cr
  24."GD.tot" kg C \cr
  25."sugar.needles" kg C \cr
  26."sugar.phloem" kg C \cr
  27."sugar.xylem.sh" kg C \cr
  28."sugar.xylem.st" kg C \cr
  29."sugar.roots" kg C \cr
  30."starch.needles" kg C \cr
  40."starch.phloem" kg C \cr
  41."starch.xylem.sh" kg C \cr
  42."starch.xylem.st"  kg C \cr
  43."starch.roots" kg C \cr)  }
}%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Documentation: Joanna Simms
%%  ~~who you are~~
}
\note{
Last updated: 27.07.2022
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 
See original "README.Rmd" file.
 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
