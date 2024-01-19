\name{Plant_N_Uptake}
\alias{Plant_N_Uptake}
\alias{plant_n_uptake}
\docType{package}
\title{
  Plant Nitrogen Uptake
}
\description{
  This function is for the soil to root uptake only.

  This function uses the nitrogen uptake functions for each of the indervidual types of N compiled for the plant.

  The demand is calculated by the decision function, then used to sclae the maximum uptake. Although demand is currently 1 for all decision functions included. (Meyer, 2010, Franklin, 2014)

  \eqn{N_to_root = (u_{organic}(FOM, T, SWC) + u_{NH4}(NH4, T, SWC) + f(NH4)u_{NO3}(NO3, T, SWC))*demand}

  The nitrogen to the root is a combination of all of the uptake functions, with the NO3 nitrogen uptake updated with a NH4 modifier. NH4, NO3 and FOM are the amounts in the soil. T is temperature in the soil B horizon and then SWC is the soil water content.

  \eqn{f(NH4) = \frac{a NH4^8}{NH4^8 + b^8}}

  Where a and b are fitted parameters and the function gives a value between 0 and 1.
}
\usage{

Plant_N_Uptake(NC_in_root_opt,
               T,
               SWC,
               m,
               NH4_in,
               NO3_in,
               FOM_in,
               N_limits_R,
               N_k_R,
               SWC_k_R,
               parameters,
               demand)

}

\arguments{
  \item{NC_in_root_opt}{
    The optimal NC ratio in the roots
  }
  \item{T}{
    Temperature
  }
  \item{SWC}{
    Soil water content
  }
  \item{m}{
    Mycorrhized potential
  }
  \item{NH4_in}{
    The NH4 in the soil
  }
  \item{NO3_in}{
    The NO3 in the soil
  }
  \item{FOM_in}{
    The nitrogen in the FOM
  }
  \item{N_limits_R}{
    Turnover of the roots that ARE mycorrhized
  }
  \item{N_k_R}{
    Parameters to control the N uptake
  }
  \item{SWC_k_R}{
    Parameters to control the N uptake
  }
  \item{parameters}{
    The parameters that control the NH4 effect on NO3 uptake
  }
  \item{demand}{
    The demand is calculated from other functions and just added here as a factor.
  }
}

\value{
    \item{N_to_plant}{The total N of all forms to the plant (N forms are not distinguished between within the plant, kg C eq}
    \item{NH4_used}{All NH4 uptaken - balance used in the symphony soil model, kg C eq}
    \item{NO3_used}{All NO3 uptaken - balance used in the symphony soil model, kg C eq}
    \item{Norg_used}{All organic N uptaken - balance used in the symphony soil model, kg C eq}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi

  Last updated 19.01.2024
}
\references{
  Meyer, Astrid, et al. "Simulating mycorrhiza contribution to forest C-and N cycling-the MYCOFON model." Plant and soil 327 2010: 493-517.

  Franklin, Oskar, et al. "Forests trapped in nitrogen limitation–an ecological market perspective on ectomycorrhizal symbiosis." New Phytologist 203.2 2014: 657-666.
}
