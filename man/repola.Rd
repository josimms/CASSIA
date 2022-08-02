\name{repola}
\alias{repola}
\title{Repola function}

\description{
  A function to impliment the repola growth equations.
}

\input{
  \item{D0}{
    Initial diameter.
  }
  \item{h0}{
    Initial height.
  }
  \item{n.year}{
    Year number in the range of years entred into the CASSIA model.
  }
  \item{needle_mas}{
    Needle mass calculated from values at Hyytiälä, method found in the CASSIA function.
  }
  \item{site}{
    Site to determine the parameters. Currently choice between "Hyde" and "Lettosuo".
  }
  \item{params}{
    The initial parameters given to the model. Default values seen in parameters_p.
  }
  \item{reps}{
    The initial values for the repola model. Default values seen in repo_p
  }
}

\output{
  repola_p, a list with four components. \cr
  1. "needle_mass" \cr
  2. "m.N.tot" \cr
  3. "m.N" \cr
  4. "m.R.tot" \cr
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
