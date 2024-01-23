\name{Symphony_Multiple_FOM_Daily}
\alias{Symphony_Multiple_FOM_Daily}
\alias{symphony_multiple_FOM_daily}
\docType{package}
\title{
  Symphony Multiple FOM Daily
}
\description{
  This model updates the soil states from CASSIA output (eventually).

  \eqn{C_{FOM}^{needles} = C_{FOM}^{needles} + L^{needles}} \cr
  \eqn{C_{FOM}^{woody} = C_{FOM}^{woody} + L^{woody}} \cr
  \eqn{C_{FOM}^{roots} = C_{FOM}^{roots} + L^{roots}} \cr
  \eqn{C_{FOM}^{mantle} = C_{FOM}^{mantle} + L^{mantle}} \cr
  \eqn{C_{FOM}^{ERM} = C_{FOM}^{ERM} + L^{ERM}}

  Where C and N are carbon and nitrogen respectively, L is letter and the compartments are refereed to explicitly. FOM is fresh organic matter and SOM is soil organic matter. These values are then aggregated for the microbe model.

  \eqn{C_{FOM} = C_{FOM}^{needles} + C_{FOM}^{woody} + C_{FOM}^{roots} + C_{FOM}^{mantle} + C_{FOM}^{ERM}} \cr
  \eqn{N_{FOM} = N_{FOM}^{needles}(N:C)^{needles} + N_{FOM}^{woody}(N:C)^{woody} + N_{FOM}^{roots}(N:C)^{roots} + N_{FOM}^{mantle}(N:C)^{mantle} + N_{FOM}^{ERM}(N:C)^{ERM}}

  where \eqn{C} and \eqn{N} are carbon and nitrogen respectively, the compartment is indicated by the superscript and \eqn{(N:C)} is the carbon nitrogen ratio of the indicated compartment.

  Then the plant and microbe usages from the previous model / iteration are calculated.

  \eqn{NH4 = NH4 - NH4_{used}^{p} - NH4_{used}^{f}} \cr
  \eqn{NO3 = NO3 - NO3_{used}^{p} - NO3_{used}^{f}} \cr
  \eqn{N^{FOM} = N^{FOM} - N_{used}^{FOM,p} - N_{used}^{FOM,f}} \cr
  \eqn{N^{SOM} = N^{SOM}}

  Where NH4_{used}^{p} is the NH4 used for the plant, f represents the fungal uptake. This is an input to the function from the mycofon model.

  TODO: carbon considerations from fungal and plant uptake!

  Next the microbe uptake is considered, and this is taken away from the N type states and C. First nitrogen balance.

  \eqn{NH4 = NH4 - u_{NH4}C^{FOM} - u_{NH4}C^{SOM}} \cr
  \eqn{NO3 = NO3 - u_{NO3}C^{FOM} - u_{NO3}C^{SOM}} \cr
  \eqn{N^{FOM} = N^{FOM} - u_{Norg}C^{FOM} - u_{Norg}C^{SOM}} \cr
  \eqn{N^{SOM} = N^{SOM} - u_{Norg}C^{SOM}}

  Where u is the uptake. Next carbon balance.

  \eqn{d_{total} = u_{C}C^{SOM} + u_{C}C^{FOM}}

  Where d is decomposition.

  TODO: The distributed \eqn{C_{FOM}} decomposition is currently based on the size of the pool rather than the type of material. This should be changed!

  \eqn{C_{FOM}^{needles} = C_{FOM}^{needles} - d_{total} \frac{C_{FOM}^{needles}}{C_{FOM}}} \cr
  \eqn{C_{FOM}^{woody} = C_{FOM}^{woody} - d_{total} \frac{C_{FOM}^{woody}}{C_{FOM}}} \cr
  \eqn{C_{FOM}^{roots} = C_{FOM}^{roots} - d_{total} \frac{C_{FOM}^{roots}}{C_{FOM}}} \cr
  \eqn{C_{FOM}^{mantle} = C_{FOM}^{mantle} - d_{total} \frac{C_{FOM}^{mantle}}{C_{FOM}}} \cr
  \eqn{C_{FOM}^{ERM} = C_{FOM}^{ERM} - d_{total} \frac{C_{FOM}^{ERM}}{C_{FOM}}}

  \eqn{C_{SOM} = C_{SOM} - u_{C}C^{microbe, SOM} + turnover*(C^{microbe, SOM} + C^{microbe, FOM})}

  \eqn{C^{microbe, FOM} = (1 + u_{SOM}C^{microbe, FOM} - r(Tmb) - turnover)C^{microbe, FOM}} \cr
  \eqn{C^{microbe, SOM} = (1 + u_{SOM}C^{microbe, SOM} - r(Tmb) - turnover)C^{microbe, SOM}}
}

\usage{

symphony_multiple_FOM_daily(Tmb,
                            SWC,
                            C_FOM_needles_old,
                            C_FOM_woody_old,
                            C_FOM_roots_old,
                            C_FOM_mantle_old,
                            C_FOM_ERM_old,
                            C_SOM_old,
                            N_SOM_old,
                            C_decompose_FOM,
                            C_decompose_SOM,
                            N_decompose_FOM,
                            N_decompose_SOM,
                            Litter_needles,
                            Litter_woody,
                            Litter_roots,
                            Litter_mantle,
                            Litter_ERM,
                            imobilisation,
                            assimilation,
                            NH4_old,
                            NO3_old,
                            NC_needles,
                            NC_woody,
                            NC_roots,
                            NC_mantle,
                            NC_ERM,
                            NH4_used_Plant,
                            NH4_used_Fungal,
                            NO3_used_Plant,
                            NO3_used_Fungal,
                            FOM_Norg_used_Plant,
                            FOM_Norg_used_Fungal,
                            SOM_Norg_used,
                            respiration_microbes_params,
                            N_limits_R,
                            N_k_R,
                            SWC_k_R,
                            NC_microbe_opt,
                            microbe_turnover)

}

\arguments{
  \item{Tmb}{
    Temperature at soil depth B
  }
  \item{SWC}{
    Soil water content
  }
  \item{C_FOM_needles_old}{
    FOM from needles
  }
  \item{C_FOM_woody_old}{
    FOM from woody
  }
  \item{C_FOM_roots_old}{
    FOM from roots
  }
  \item{C_FOM_mantle_old}{
    FOM from mycelium
  }
  \item{C_FOM_ERM_old}{
    FOM from ERM
  }
  \item{C_SOM_old}{
    SOM
  }
  \item{N_SOM_old}{
    N in SOM
  }
  \item{C_decompose_FOM}{
    C in the FOM decomposers
  }
  \item{C_decompose_SOM}{
    C in the FOM decomposers
  }
  \item{N_decompose_FOM}{
    N in the FOM decomposers
  }
  \item{N_decompose_SOM}{
    N in the FOM decomposers
  }
  \item{Litter_needles}{
    The input of needle litter from CASSIA
  }-
  \item{Litter_woody}{
    The input of woody litter from CASSIA
  }
  \item{Litter_roots}{
    The input of root litter from CASSIA
  }
  \item{Litter_mantle}{
    The input of mantle litter from CASSIA
  }
  \item{Litter_ERM}{
    The input of ERM litter from CASSIA
  }
  \item{NH4_old}{
    NH4 from the soil
  }
  \item{NO3_old}{
    NO3 from the soil
  }
  \item{NC_needles}{
    The nitrogen : carbon ratio for the needle FOM.
  }
  \item{NC_woody}{
    The nitrogen : carbon ratio for the needle FOM.
  }
  \item{NC_roots}{
    The nitrogen : carbon ratio for the needle FOM.
  }
  \item{NC_mantle}{
    The nitrogen : carbon ratio for the needle FOM.
  }
  \item{NC_ERM}{
    The nitrogen : carbon ratio for the needle FOM.
  }
  \item{NH4_used_Plant}{
    NH4 used by / allocated by the plant
  }
  \item{NH4_used_Fungal}{
    NH4 used by / allocated by the fungus
  }
  \item{NO3_used_Plant}{
    NO3 used by / allocated by the plant
  }
  \item{NO3_used_Fungal}{
    NO3 used by / allocated by the fungus
  }
  \item{FOM_Norg_used_Plant}{
    NO3 used by / allocated by the plant
  }
  \item{FOM_Norg_used_Fungal}{
    NO3 used by / allocated by the fungus
  }
  \item{SOM_Norg_used}{
    SOM used by / allocated by the plant
  }
  \item{respiration_microbes_params}{
    Parameters that control the microbial respiration
  }
  \item{N_limits_R}{
    The parameters that change the uptake parameters for microbes
  }
  \item{N_k_R}{
    The parameters that change the uptake parameters for microbes
  }
  \item{SWC_k_R}{
    The parameters that change the uptake parameters for microbes
  }
  \item{NC_microbe_opt}{
    Optimal N:C ratio in the microbes
  }
  \item{microbe_turnover}{
    The turnover of microbes
  }
}

\value{
  \item{C_decompose_FOM}{
    C in the FOM decomposers
  }
  \item{C_decompose_SOM}{
    C in the SOM decomposers
  }
  \item{C_FOM_ERM}{
    C in the FOM ERM biomass.
  }
  \item{C_FOM_mantle}{
    C in the FOM mantle biomass.
  }
  \item{C_FOM_needles}{
    C in the FOM needle biomass.
  }
  \item{C_FOM_roots}{
    C in the FOM roots biomass.
  }
  \item{C_FOM_woody}{
    C in the FOM woody biomass.
  }
  \item{C_SOM}{
    C in the SOM biomass.
  }
  \item{N_decompose_FOM}{
    N in the decompoers of FOM
  }
  \item{N_decompose_SOM}{
    N in the decompoers of SOM
  }
  \item{F_SOM}{
    N in the FOM biomass.
  }
  \item{N_SOM}{
    N in the SOM biomass.
  }
  \item{NC_ERM}{
    NC Ratio in the ERM
  }
  \item{NC_mantle}{
    NC Ratio in the mantle
  }
  \item{NC_needles}{
    NC Ratio in the needles
  }
  \item{NC_roots}{
    NC Ratio in the roots
  }
  \item{NC_woody}{
    NC Ratio in the woody biomass
  }
  \item{NH4}{
    NH4 in soil.
  }
  \item{NO3}{
    NO3 in soil
  }
  \item{SOM_Norg_used}{
    The Norg in SOM decomposed to get the extra nitrogen that the microbes want.
  }
  \item{Microbe_respiration}{
    The unit mass respiration for microbes.
  }
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi

  Last updated 20.01.2024
}
\references{
  Perveen, Nazia, et al. "Priming effect and microbial diversity in ecosystem functioning and response to global change: a modeling approach using the SYMPHONY model." Global change biology 20.4 (2014): 1174-1190.
}
