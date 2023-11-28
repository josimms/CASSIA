\name{Mycofon_Balence}
\alias{mycofon_balence}
\alias{mycofon_balence}
\docType{package}
\title{
  Mycofon model (todo: reference), with updates to include the uptake properties for three types of nitrogen as well as the decision functionality.
}
\description{
  This function brings together the C or N balance as in the original MYCOFON model, the decision functions for resource transfer, fungal growth and multiple N uptake functions. The equations in this function are below.
  
  \eqn{m = \frac{C_f}{B^{rf}_{opt}C^r}}
  
  Where N is nitrogen, C is carbon, r roots, f fungai, B biomass, opt is optimum.
  
  The uptake and decision functions are then calulated from the balance of the last iteration. Note that the plant decision requires the N alllocated, so the mycorrhizal decision is calculated first.
  
  Then the balance is updated again.
  
  \eqn{C^r = C^r + C^{CASSIA}_{Allocated} - (1 - m)t^rC^r - mt^mC^r - r(Tmb)C^r - d^r(C^{CASSIA})}

  Where, C^{CASSIA}_{Allocated} is the carbon allocated from the CASSIA model, t is turnover, which is different for micorhized (m) and non-mycorrhized roots (r), Tmb is temperature of the soil, d() is decision and r() is respiration. Growth is currently handled in CASSIA, and therefore has not been added here!
  
  \eqn{C^f = C^f + d^r(C^{CASSIA}) - 0.5 t^{mantle} C^{mantle} - 0.5 t^{ERM} C^{ERM} - r(Tmb)  C^{f}}
  
  \eqn{N^r = N^r + d^f(N^{f}) + u^r(T, SWC, NH4, NO3, FOM_Norg)C^r - (1 - m) t^r R^r_{C:N} N^r - mt^m R^r_{C:N} N^r - R^r_{C:N} N^{CASSIA}_{Allocated}}
  
  R^r_{C:N} fungal NC ratio
  
  \eqn{N^f = N^f + u^f(T, SWC, NH4, NO3, FOM_Norg)C^f R^f_{C:N} - t^{ERM} C^{ERM} R^{ERM}_{C:N} - t^{mantle} C^{mantle} R^{mantle}_{C:N} - d^f(N^f)}
}
\usage{

mycofon_balence(C_roots,
                N_roots,
                optimal_root_fungal_biomass_ratio,
                C_fungal,
                N_fungal,
                turnover_roots,
                turnover_roots_mycorrhized,
                turnover_mantle,
                turnover_ERM,
                respiration_parameters_R,
                NH4,
                NO3,
                FOM_Norg,
                NC_in_fungai_opt,
                T,
                Tsb,
                SWC,
                N_limits_Plant,
                N_k_Plant,
                SWC_k_Plant,
                N_limits_Fungal,
                N_k_Fungal,
                SWC_k_Fungal,
                NC_fungal_opt,
                mantle_mass,
                ERM_mass,
                parameters_NH4_on_NO3,
                growth_C,
                growth_N,
                max_C_allocation_CASSIA,
                allocation_N_to_rest_of_plant,
                mycofon_stratergy)
}

\arguments{
  \item{C_roots}{
    Carbon in roots.
  }
  \item{N_roots}{
    Nitrogen in roots.
  }
  \item{optimal_root_fungal_biomass_ratio}{
    The ratio of root and fungal biomass, as in mochrizzal rate.
  }
  \item{C_fungal}{
    Carbon in mycorrhzia.
  }
  \item{N_fungal}{
    Nitrogen in mycorrhiza.
  }
  \item{turnover_roots}{
    Turnover of the roots that ARE NOT mycorrhized
  }
  \item{turnover_mycorrhized}{
    Turnover of the roots that ARE mycorrhized
  }
  \item{turnover_mantle}{
    Turnover of the mantle
  }
  \item{turnover_ERM}{
    Turnover of the ERM
  }
  \item{respiration_parameters_R}{
    Respiration parameters
  }
  \item{NH4}{
    Amount of NH4 in the soil
  }
  \item{NO3}{
    Amount of NO3 in the soil
  }
  \item{FOM_Norg}{
    The nitrogen in FOM
  }
  \item{NC_in_fungai_opt}{
    The idea NC ratio in the fungus
  }
  \item{T}{
    Temperature of the air, 'C
  }
  \item{Tsb}{
    Temperature of the B horizon of the soil, 'C
  }
  \item{SWC}{
    Soil water content, \%
  }
  \item{N_limits_Plant}{
    Parameters to control the N uptake
  }
  \item{N_k_Plant}{
    Parameters to control the N uptake
  }
  \item{SWC_k_Plant}{
    Parameters to control the N uptake
  }
  \item{N_limits_Fungal}{
    Parameters to control the N uptake
  }
  \item{N_k_Fungal}{
    Parameters to control the N uptake
  }
  \item{SWX_k_Fungal}{
    Parameters to control the N uptake
  }
  \item{mantle_mass}{
    The mass of the mantle section of the fungal mass
  }
  \item{ERM_mass}{
    The mass of the ERM section of the fungal mass 
  }
  \item{parameters_NH4_on_NH3}{
    The additional parameters to control the uptake of the NO3 from NH4 in the plant uptake.
  }
  \item{growth_C}{
    The growth parameter for the C limitation on fungal growth.
  }
  \item{growth_N}{
    The growth parameter for the N limitation on fungal growth.
  }
  \item{max_C_allocation_CASSIA}{
    The maximum allocation of C from CASSIA to mycorrhiza
  }
  \item{allocation_N_to_rest_of_plant}{
    The allocation of the root N to the rest of the plant. In this case to the CASSIA model. NOTE, when the CASSIA links are set better, then this won't be an input.
  }
  \item{mycofon_stratergy}{
    TURE or FALSE value, Franklin is the default stratergy.
  }
}

\value{
  \item{C_roots}{Carbon in the roots after the balance, kg C eq.}
  \item{C_fungal}{Carbon in the fungi after the balance, kg C eq.}
  \item{N_roots}{Nitrogen in the roots after the balance, kg C eq.}
  \item{N_fungal}{Nitrogen in the fungi after the balance, kg C eq.}
  \item{uptake_plant}{The amount of N the plant will uptake}
  \item{uptake_fungal}{The amount of N the fungi will uptake}
  \item{from_CASSIA}{The amount of C given from CASSIA})
  \item{to_CASSIA}{The amount of N given to the rest of the plant i.e. CASSIA}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi
}
\references{
  Meyer, Astrid, et al. "Simulating mycorrhiza contribution to forest C-and N cycling-the MYCOFON model." Plant and soil 327 2010: 493-517.
  
  Franklin, Oskar, et al. "Forests trapped in nitrogen limitation–an ecological market perspective on ectomycorrhizal symbiosis." New Phytologist 203.2 2014: 657-666.
}
