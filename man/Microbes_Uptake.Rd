\name{Microbe_Uptake}
\alias{Microbe_Uptake}
\alias{micorbe_uptake}
\docType{package}
\title{
  Micorbe Uptake
}
\description{
  This uptake is taken from the Symphony model, considering the N uptake. The original C limitation is changed to a organic nitrogen limitation, assuming that the uptake of organic nitrogen is limiting C decomposition.

  The uptake follows the form.

  \eqn{ (a r(T), C_{microbe} + b - a)  \min \{\text{organic nitrogen limitation}, \text{nitrogen limitation} \}}

  Where a is the microbe N:C and b is the soil compartment N:C ratio.

  Organic nitrogen limitation = \eqn{u(Norg)} \cr
  Nitrogen limitation = \eqn{\frac{i(u(NH4) + u(NO3)) + (N:C)^{microbes}_{opt} r(T) C^{microbe}}{b - a}}

  Where u() is the uptake function, \eqn{(N:C)^{microbe}_{opt}} is the N:C ratio and r(T) is the microbe respiration.

  The outputs are then scaled by the lowest limitation.

  uptake * \eqn{\left( \frac{\text{min limitation}}{\text{max limitation}} \right)}
}

\usage{
  Microbe_Uptake(C_microbe,
                 N_micorbe,
                 NC_microbe_opt,
                 NH4_avaliable,
                 NO3_avaliable,
                 Norg_avaliable,
                 T,
                 SWC,
                 N_limits_R,
                 N_k_R,
                 SWC_k_R,
                 SOM_decomposers,
                 Norg_avaliable_FOM)
}

\arguments{
  \item{C_microbe}{
    The C in a microbe
  }
  \item{N_microbe}{
    The N in a microbe
  }
  \item{C_soil_compartment}{
    C in the relevent soil compartment
  }
  \item{NC_microbe_opt}{
    The N:C microbe ratio
  }
  \item{NH4_avaliable}{
    The NH4 in the soil
  }
  \item{NO3_avalibale}{
    The NO3 in the soil
  }
  \item{Norg_avaliable}{
    The Norg in the FOM
  }
  \item{T}{
    Tempature at the soil horizon
  }
  \item{SWC}{
    Soil water content
  }
  \item{NC_Litter}{
    The N:C litter ratio
  }
  \item{a}{
    Assimilation
  }
  \item{i}{
    Imoblisation
  }
  \item{assimilation}{
    Decopmosition of SOM by SOM decomposers.
  }
  \item{N_limits_R}{
    Parameters to control the N uptake
  }
  \item{N_k_R}{
    Parameters to control the N uptake
  }
  \item{SWC_k_R}{
    Parameters to control the N uptake
  }
  \item{SOM_decomposers}{
    Boolian values if it is a SOM decomposer that should be calculated
  }
  \item{respiration_microbes_params}{
    Respiration parameters
  }
}

\value{
    \item{NH4_uptaken}{NH4 uptaken by relevant microbes (kg C eq)}
    \item{NO3_uptaken}{NO3 uptaken by relevant microbes (kg C eq)}
    \item{Norg_uptaken}{Norg uptaken by relevant microbes (kg C eq)}
    \item{C_uptaken}{C uptaken by relevant microbes (kg C eq)}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi
}
\references{
  Perveen, Nazia, et al. "Priming effect and microbial diversity in ecosystem functioning and response to global change: a modeling approach using the SYMPHONY model." Global change biology 20.4 (2014): 1174-1190.
}

