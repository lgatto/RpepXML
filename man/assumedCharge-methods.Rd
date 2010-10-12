\name{assumedCharge-methods}
\docType{methods}
\alias{assumedCharge-methods}
\alias{assumedCharge<-}
\alias{assumedCharge}

\title{ Accessor methods of the \code{assumedCharge} slot }

\description{
  These methods allow to retrieve or update the \code{assumedCharge} slot
  of an instance of \code{"\linkS4class{SearchResult}"}. The
  \code{assumedCharge} of all the search results of an
  \code{"\linkS4class{MSMSpepXML}"} can also be retrieved. The values
  are set automatically when creating an instance of
  \code{"\linkS4class{MSMSpepXML}"} with \code{\link{parseMSMSpepXML}},
  and should generally not need to be changed.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SpectrumQuery")}}{
      object of class \code{"\linkS4class{SpectrumQuery}"} is given as
	parameter and a numeric value is returned. }
    \item{\code{signature(object = "MSMSpepXML")}}{
      object of class \code{"\linkS4class{MSMSpepXML}"} is given as
      parameter and a numeric of length equal to the number of spectrum
      queries is returned.
    }

  }
}

\keyword{methods}
