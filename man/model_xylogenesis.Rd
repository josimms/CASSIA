\name{model_xylogenesis}
\alias{model_xylogenesis}
\title{Xylogenesis model function}

\description{
  Function that contains the xylogenesis effect of the model. \cr

  Not yet integrated into the CASSIA model. \cr

  TODO: should have a deeper explaination when it is actualy included in the model and when the function itself is finished.
}
\arguments{
 \item{g.sD.T}{
   ??
 }
  \item{LD}{
   Daily secondary growth
 }
  \item{fD}{
   Reduction function for the daily growth.
 }
  \item{storage_reduction}{
   TRUE or FALSE . The reduction applied if there is a storage reduction.
 }
 \item{environment_effect_xylogenesis}{
   TRUE or FALSE . Determines growth considering the environmental effect or not.
   If FALSE then...
 }
 \item{Temp}{
   Temperature 'C
 }
}

\value{
  Output is in the form of a list.
  \item{GD}{
   Daily potential number of new cells per dat (in one radial cell row)
  }
  \item{tot.cells}{
   Cumulative sum of the daily potential number of new cells.
  }
  \item{n.E}{
   Number of cells in enlargement.
  }
  \item{n.W}{
    Number of cells in wall formation.
  }
  \item{n.M}{
    Number of cells in the mature phase.
  }
  \item{ew_cells}{
    Early wood cells.
  }
  \item{lw_cells}{
    Late wood cells.
  }
  \item{ring_width}{
    The ring width, daily time series, mm
  }
  \item{ew_width}{
    The early wood width, daily time series, mm
  }
  \item{lw_width}{
    The late wood width, daily time series, mm
  }
  \item{en.growth}{
    Carbon to enlargement per day, kg C day-1
  }
  \item{en.release}{
    The carbon used for enlargement is release after a delay cause by tau.Ee, kg C day-1
  }
  \item{wall.growth}{
    The use of carbon in wall growth, kg C day-1
  }
  \item{wall.tot}{
    Cumulative use of carbon in wall growth. kg C year-1
  }
  \item{cell_d_final}{
    Daily diameter, daily values. Formed by repetition of the cell diameter for early wood and late wood, for calculated number of cells.
  }
  \item{cell_wall_thickness}{
    Cell wall thickness, daily values. Formed by repetition of the wall thickness for early wood and wall tickeness for late wood, for calculated number of cells.
  }
  \item{cell_density}{
    Cell wall density, daily values. Formed by repetition of the wall density for early wood and wall density for late wood, for calculated number of cells.
  }
  \item{ring_density}{
    Carbon density of the rings by volume. Worked out with the cumulative sum of CW at the end point divided by the last number of cells multiplied by their diameter squared.
  }
}

\references{

}

\author{
  Model: Pauliina Schiestl-Aalto.
  Documentation: Joanna Simms.
}

\note{
  Created: 27.07.2022
  Last updated: 27.07.2022
}
