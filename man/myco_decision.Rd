\name{myco_decision}
\alias{Myco Decision}
\docType{package}
\title{
  Myco Decision
}
\description{
  Currently two strategies are used. These are derived from Fraklin et al 2014 and the Mycofon model (Meyer, 2010).

  The stratergy of Mycofon gives;

  a demand of 1

  and

  Transfer = \eqn{\max \{N_f(1 - \frac{\frac{N_r}{C_r}}{(N:C)^{*}}), 0.0 \}}

  where \eqn{N_f}, \eqn{N_r} are the N in the fungi and root respectively, \eqn{C_r} is the C in the root and (N:C)^{*} is the optimal N:C ration in the root.

  The strategy of Franklin gives;

  a demand of 1

  and

  Transfer = \eqn{\max \{N_f - G_f, 0.0 \}}

  Where N and G are the nitrogen and growth of the fungi.

}

\usage{
myco_decision(C_fungal,
              N_fungal,
              C_roots,
              N_roots,
              NC_fungal_opt,
              growth_C,
              growth_N)
}

\arguments{
  \item{C_fungal}{
    The C in the fungi.
  }
  \item{N_fungal}{
    The N in the fungi.
  }
  \item{C_roots}{
    The C in the roots.
  }
  \item{N_roots}{
    The N in the roots.
  }
  \item{NC_fungal_opt}{
    The optimal N:C ratio in the fungi.
  }
  \item{growth_C}{
    The fungal growth C component factor.
  }
  \item{growth_N}{
    The fungal growth N component factor.
  }
}

\value{
  \item{Mycofon_demand}{1}
  \item{Mycofon_allocarion}{The N allocated to the plant by mycorrhiza according to the mycofon method}
  \item{Franklin_demand}{1}
  \item{Franklin_allocarion}{The N allocated to the plant by mycorrhiza according to the Franklin method}
}

\author{
  Joanna Simms, joanna.x.simms@helsinki.fi
}
\references{
  Meyer, Astrid, et al. "Simulating mycorrhiza contribution to forest C-and N cycling-the MYCOFON model." Plant and soil 327 2010: 493-517.

  Franklin, Oskar, et al. "Forests trapped in nitrogen limitation–an ecological market perspective on ectomycorrhizal symbiosis." New Phytologist 203.2 2014: 657-666.
}
