\name{Hyytiala_Data_Creation}
\alias{Hyytiala_Data_Creation}

\title{Hyytiala_Data_Creation function. A function which downloads the data from the Hyytiala, SMEAR II database}

\arguments{
  \item{raw.directory}{
    This need s directory to download the files from the internet and to create a file for the weather data, after the data has been created in the format needed for CASSIA.
  }
  \item{year_start}{
    The year that the data download starts from
  }
  \item{year_end}{
    The year that the data download ends
  }
  \item{download}{
    If you don't want the download section (as the files are already downloaded), you can toggle this T or F.
  }
  \item{clean_data}{
    If you don't want to clean the data (for example, if you only want to use this to download the data), you can toggle this T or F.
  }
}

\value{
  Outputs a dataframe with the format / headers needed for CASSIA. \cr
  Note! You can also look at the code to see the formatting needed for CASSIA to work.
}

\references{
  Data from: https://smear.avaa.csc.fi/
}

\author{
  Documentation / Code: Joanna Simms
}

\note{
  Created: 05.05.2023
  Last updated: 05.05.2023
}
