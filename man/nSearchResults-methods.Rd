\name{nSearchResults-methods}
\docType{methods}
\alias{nSearchResults}
\alias{nSearchResults-methods}

\title{ Returns the number of search results. }

\description{
  These methods return the number of search results that have been 
  identified. The class of the returned object will depend on the
  argument provided, as described below. Note that simplification is
  attempted by 'sapply'. 
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "SpectrumQuery")}}{
      an object of class \code{"\linkS4class{SpectrumQuery}"} is used as
      parameter and an integer is returned.
    }
    \item{\code{signature(object = "MSMSpepXML")}}{
      an object of class \code{"\linkS4class{MSMSpepXML}"} is used as
      parameter and an integer of length equal to the number of spectrum
      queries is returned.
    }
  }
}

\keyword{methods}

