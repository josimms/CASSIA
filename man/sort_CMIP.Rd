\name{sort_CMIP}
\alias{sort_CMIP}

\title{Sort CMIP data.}

\description{
  Data prepared for the nitrogen phydro photosynthesis model.
  
  The data should be predownloaded.
  
  Function code is found in CMIP6.R
}

\arguments{
  \item{TEMP_DIR}{
    The place the CMIP6 data is stored.
  }
  \item{ANSWER_DIR}{
    The path the sorted data is saved.
  }
  \item{path_test}{
    Data for the bias correction.
  }
}

\value{
  There isn't a return for this function. The monthly and daily datasets should be found in ANSWER_DIR. 
}

\references{
  Data from: https://pcmdi.llnl.gov/CMIP6/
}

\author{
  Documentation / Code: Joanna Simms
}

\note{
  Created: 15.10.2024
  Last updated: 15.10.2024
}
