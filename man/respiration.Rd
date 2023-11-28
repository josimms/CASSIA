\name{respiration}
\alias{Respiration}
\docType{package}
\title{
  Respiration
}
\description{
  The respiration is the same for all of the organisms, but the parameters are different.
  
  \eqn{\max \{Rm \exp(\frac{Q10 T}{10}) - \exp(- \frac{1}{2} \log(Q10)), 0.0 \}}
  Where Rm and Q10 are photosynthesis parameters. T is temperature. The formula is from the CASSIA formation (Schistel-Aalto, 2015)
}

\usage{
respiration(T,
            a,
            b)
}

\arguments{
  \item{T}{
    Temperature of the air or the soil depending on the microbes.
  }
  \item{a}{
    Scaling factor.
  }
  \item{b}{
    Exponential factor.
  }
}

\value{
  \item{output}{Amount of the respiration, kg C eq}
}

\author{
  Joanna Simms, joanna.x.simms@ad.helsinki.fi
}
\references{
  Schiestl‐Aalto, Pauliina, et al. "CASSIA–a dynamic model for predicting intra‐annual sink demand and interannual growth variation in S cots pine." New Phytologist 206.2 (2015): 647-659.
}
