\name{Myco_Growth}
\alias{myco_growth}
\docType{package}
\title{
  Mycorrhizal Growth
}
\description{
  The growth of mycorrhiza is controlled by the limitation in C or N.
  
  \eqn{\min \{a C, b N \}}
  
  Where a is the growth rate for fungi dependent on C, and b is the growth rate of fungi dependent on N.
  
  Linear functions as in Meyer 2010 and Franklin 2014.
}

\usage{
respiration(C_fungal,
            N_fungal,
            a,
            b)
}

\arguments{
  \item{C_fungal}{
    Carbon in the fungal biomass
  }
  \item{N_fungal}{
    Nitrogen in the fungal biomass
  }
  \item{a}{
    Factor for the C decomposition
  }
  \item{b}{
    Factor for the N decomposition
  }
}

\value{
  \item{growth}{Mycorrhizal growth in terms of kg C.}
  \item{carbon_used}{Cabon used for growth in terms of kg C.}
  \item{nitrogen_used}{Nitrogen used for growth in terms of kg C eq.}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi
}
\references{
  Meyer, Astrid, et al. "Simulating mycorrhiza contribution to forest C-and N cycling-the MYCOFON model." Plant and soil 327 2010: 493-517.
  
  Franklin, Oskar, et al. "Forests trapped in nitrogen limitation–an ecological market perspective on ectomycorrhizal symbiosis." New Phytologist 203.2 2014: 657-666.
}
