\name{precNeutralMass-methods}
\docType{methods}
\alias{precNeutralMass-methods}
\alias{precNeutralMass<-}
\alias{precNeutralMass}

\title{ Returns the precursor's neutral mass. }

\description{
  This method allows to retrieve and set the neutral mass of the
  precursor ion used for MSMS peptide identification. The values are set
  automatically when creating an instance of
  \code{"\linkS4class{MSMSpepXML}"} with \code{\link{parseMSMSpepXML}},
  and should generally not need to be changed. The method can also be
  used to retrieve the precursor's neutral masses of all the spectrum
  queries of an \code{"\linkS4class{MSMSpepXML}"} object.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SpectrumQuery")}}{
      an object of class \code{"\linkS4class{SpectrumQuery}"} is used as
      parameter and an numeric is returned.
    }
    \item{\code{signature(object = "MSMSpepXML")}}{
      an object of class \code{"\linkS4class{MSMSpepXML}"} is used as
      parameter and an numeric or length equal to the number of spectrum
      queries is returned.    
    }
  }
}
\keyword{methods}

