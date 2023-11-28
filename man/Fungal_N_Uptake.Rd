\name{Fungal_N_Uptake}
\alias{Fungal_N_Uptake}
\alias{fungal_n_uptake}
\docType{package}
\title{
  Fungal Nitrogen Uptake
}
\description{
  This function uses the nitrogen uptake functions for each of the indervidual types of N compiled for the fungae.
  
  The demand is calculated by the decision function, although is currently 1 for all decision functions included. (Meyer, 2010, Franklin, 2014)
  
  Then the output is the N uptake of each of the types of nitrogen calculated in parallel (Meyer, 2010) a then downscaled with the demand function (me).
  
  \eqn{u_{total} = (u_{Org} + u_{NO3} + u_{NH4})*demand}.
  
  Where \eqn{u_{total}} is the total fungal uptake. \eqn{u_{Org}}, \eqn{u_{NO3}} and \eqn{u_{NH4}}, are the fungal uptake of each of the forms of nitrogen. And demand is an input to the function caluclated in the myco_decision function.
  
}

\usage{

Fungal_N_Uptake(C_fungal,
                N_fungal,
                NC_fungal_opt,
                mantle_mass,
                ERM_mass,
                percentage_C_biomass,
                T,
                SWC,
                NH4,
                NO3,
                FOM_Norg,
                N_limits_R,
                N_k_R,
                SWC_k_R)

}

\arguments{
  \item{C_fungal}{
    C in the fungal mass
  }
  \item{N_fungal}{
    N in the fungal mass
  }
  \item{NC_fungal_opt}{
    The optimal ratio between C and N 
  }
  \item{mantle_mass}{
    The fungal mass that is mantle
  }
  \item{ERM_mass}{
    The fungal mass that is ERM
  }
  \item{percentage_C_biomass}{
    The percentage of carbon in biomass
  }
  \item{T}{
    The temperature at the soil horizon B
  }
  \item{SWC}{
    The soil water content
  }
  \item{NH4}{
    NH4 in the soil
  }
  \item{NO3}{
    NO3 in the soil
  }
  \item{FOM_Norg}{
    The nitrogen in the FOM
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
}

\value{
    \item{N_to_fungal}{The total N of all forms to the fungi (N forms are not distinguished between within the fungi, kg C eq}
    \item{NH4_used}{All NH4 uptaken - balance used in the symphony soil model, kg C eq}
    \item{NO3_used}{All NO3 uptaken - balance used in the symphony soil model, kg C eq}
    \item{Norg_used}{All organic N uptaken - balance used in the symphony soil model, kg C eq}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi
}
\references{
    Meyer, Astrid, et al. "Simulating mycorrhiza contribution to forest C-and N cycling-the MYCOFON model." Plant and soil 327 2010: 493-517.
  
  Franklin, Oskar, et al. "Forests trapped in nitrogen limitation–an ecological market perspective on ectomycorrhizal symbiosis." New Phytologist 203.2 2014: 657-666.
}
