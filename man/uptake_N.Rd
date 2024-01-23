\name{uptake_N}
\alias{Uptake_N}
\alias{uptake_n}
\docType{package}

\title{
  Uptake Organic
}
\description{
  The organic uptake of N is currently controlled by considering the maximum possible uptake and then this is downsclaed with functions, taking a value between 0 and 1. These functions take the same logic, but not the same form as the PRELES equations.

  Basic uptake

  \eqn{u(N) = \frac{k N^8}{N_{limit}^8 + N^8}}
  Where N is the concentration of a form of nitrogen, N_limit and k are parameters to control the uptake.

  Temperature effect

  \eqn{f_{T}(T) = \frac{T+20}{55}}
  Where T is the temperature, 20 and 55 are parameters to control the uptake.

  Water effect

  \eqn{f_{SWC}(SWC) = \frac{k SWC^8}{0.5^8 + SWC^8}}
  Where SWC is soil water content and 0.5 and k are parameters.

  This comes together to make the basic output:

  \eqn{u_{actual} = f_T(T)f_{SWC}(SWC)u(N)}.
}

\usage{
uptake_n(N,
        T,
        N_limit,
        k,
        SWC,
        SWC_sat)
}

\arguments{
  \item{N}{
    Nitrogen, \eqn{\mu} M
  }
  \item{T}{
    Temperature (usually the temperature at soil horizon B is used). degrees C
  }
  \item{N_limit}{
    Factor for the pure N uptake. \eqn{\mu} M
  }
  \item{k}{
    Factor for the pure N uptake.
  }
  \item{SWC}{
    Soil water concentration, \%
  }
  \item{SWC_sat}{
    Saturation point for the SWC effect on N uptake. \%
  }
}

\value{
    \item{N uptake}{The uptake of a form of N, \eqn{\mu} M}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi

  Last updated 19.01.2024
}
\references{
  Näsholm, Torgny, Knut Kielland, and Ulrika Ganeteg. "Uptake of organic nitrogen by plants." New phytologist 182.1 (2009): 31-48.

  Data fitting for automatic parameters;
  Oyewole, Olusegun Ayodeji, et al. "Soil diffusive fluxes constitute the bottleneck to tree nitrogen nutrition in a Scots pine forest." Plant and Soil 399 (2016): 109-120.
}
