\name{nSearchHits-methods}
\docType{methods}
\alias{nSearchHits}
\alias{nSearchHits-methods}

\title{ Returns the number of search hits. }
\description{
  These methods return the number of search hits that have been
  identified for one or several search results. The class of the
  returned object will depend on the argument provided, as described
  below. Note that simplification is attempted by 'sapply'.
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchResult")}}{
      an object of class \code{"\linkS4class{SearchResult}"} is used as
      parameter and an integer is returned. 
    }
    \item{\code{signature(object = "SpectrumQuery")}}{
      an object of class \code{"\linkS4class{SpectrumQuery}"} is used as
      parameter and an integer of length equal to the number of search
      results is returned.
    }
    \item{\code{signature(object = "MSMSpepXML")}}{
      an object of class \code{"\linkS4class{MSMSpepXML}"} is used as
      parameter and a list is returned.

    }
  }
}

\keyword{methods}

