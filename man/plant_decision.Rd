\name{plant_decision}
\alias{Plant Decision}
\docType{package}
\title{
  Plant Decision
}
\description{
  Currently two strategies are used. These are derived from Fraklin et al 2014 and the Mycofon model (Meyer, 2010).

  The strategy of Mycofon gives;

  a demand of 1

  and

  transfer = \eqn{\max \{allo_{max} C_r, allo(C_r - (R:F)*) - C_f \}}

  where \eqn{allo_{max} = 1 - (1 - [exp(-50*N_r)]^3)} if \eqn{N_r < 0.01}, else it's 0.2 and

  \eqn{allo = N_{allo} / (N_r + N_{allo})} when \eqn{N_{allo} < 0.5 N_r} else it's 1.

  \eqn{C_r} and \eqn{C_f} are the carbon in the roots and fungi respectively, \eqn{(R:F)*} is the optimal root fungal biomass ratio, \eqn{N_r} is the nitrogen in the roots and N_{allo} is the nitrogen allocated by the mycorrhiza.

  The strategy of Franklin gives;

  a demand of 1

  and

  transfer = \eqn{\max \{ \max(CASSIA, 0.0) \}}

  where CASSIA is the maximum C allocation to the mycorrhiza calculated in CASSIA.
}

\usage{
plant_decision(C_roots,
               N_roots,
               C_fungal,
               optimal_root_funga_biomass_ratio,
               N_allo,
               max_C_allocation_CASSIA)
}

\arguments{
  \item{C_roots}{
    The C in the roots.
  }
  \item{N_roots}{
    The N in the roots.
  }
  \item{C_fungi}{
    The C in the fungi.
  }
  \item{optimal_root_funga_biomass_ratio}{
    The optimum value of the root:fungai biomass ratio.
  }
  \item{N_allo}{
    N allocated to the plant
  }
  \item{max_C_allocation_CASSIA}{
    Maximum C allocated to the mycorrhiza from CASSIA
  }
}

\value{
  \item{Mycofon_demand}{1}
  \item{Mycofon_allocarion}{The C allocated to the fungi by mycorrhiza according to the mycofon method}
  \item{Franklin_demand}{1}
  \item{Franklin_allocarion}{The C allocated to the fungi by mycorrhiza according to the Franklin method}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi
}
\references{
  Meyer, Astrid, et al. "Simulating mycorrhiza contribution to forest C-and N cycling-the MYCOFON model." Plant and soil 327 2010: 493-517.

  Franklin, Oskar, et al. "Forests trapped in nitrogen limitation–an ecological market perspective on ectomycorrhizal symbiosis." New Phytologist 203.2 2014: 657-666.

}
