\name{PRELES_GPP}
\alias{PRELES_GPP}

\title{PRELES_GPP function. An adaptation of the PRELES function from the Rprebasso package}

\arguments{
  \item{photoparameters}{
    Options of p values for the PRELES function and the growth_photo_coef. \cr
    1 uses the normal vector, growth_photo_coef = 1 \cr
    2 uses c(rep(NA,4), 0.8957, rep(NA, 4), 0.03473, rep(NA, 20)), growth_photo_coef = 1.24 \cr
    3 uses c(rep(NA,4), 0.9357, rep(NA, 4), 0.03273, rep(NA, 20)), growth_photo_coef = 1.35 \cr
    3 uses c(rep(NA,4), 0.9357, rep(NA, 4), 0.03273, rep(NA, 20)), growth_photo_coef = 1 \cr
  }
  \item{Temp}{
    Daily mean temperature. degree C
  }
  \item{PF}{
    Daily photosynthesis from a model such as SPP. kg C
  }
  \item{Tsa}{
    Daily mean temperature for the soil at horizon A. degree C
  }
  \item{Tsb}{
    Daily mean temperature for the soil at horizon B. degree C
  }
  \item{M.soil}{
    Daily soil moisture. m3 m-3
  }
  \item{Rain}{
    Daily sum of precipitation. mm
  }
}

\value{
  The outputs are a list of daily photosynthesis calculated by PRELES with the determined conditions and the growth_photo_coefficient.
}

\references{
  # TODO: add a reference
}

\author{
  Documentation: Joanna Simms
  Content: Pauliina Schiestl-Aalto
}

\note{
  Created: 27.07.2022
  Last updated: 27.07.2022
}
