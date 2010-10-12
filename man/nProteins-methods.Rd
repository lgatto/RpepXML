\name{nProteins-methods}
\docType{methods}
\alias{nProteins}
\alias{nProteins-methods}

\title{ Returns the number of identified proteins. }

\description{
 These methods return the number of proteins that have been identified
 for any individual seach hit. The class of the returned object will
 depend on the argument provided, as described below. Note that
 simplification is attemplted by 'sapply'.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchHit")}}{
      an object of class \code{"\linkS4class{SearchHit}"} is used as
      parameter and a single integer value is returned.
    }
    \item{\code{signature(object = "SearchResult")}}{
      an object of class \code{"\linkS4class{SearchResult}"} is used as
      parameter and an integer is returned of length equal to the number
      of search hits in the search result. If the search result has no
      search hits, 0 is returned. 
    }
    \item{\code{signature(object = "SpectrumQuery")}}{
      an object of class \code{"\linkS4class{SpectrumQuery}"} is used as
      parameter the output will depend on the structure of the data. A
      spectrum without search results will return 0. If one search
      result exists, the returned value will be an integer or a matrix,
      depending on the number of search hits. If the sectrum query
      contains several search results, a list is returned.
    }
    \item{\code{signature(object = "MSMSpepXML")}}{
      an object of class \code{"\linkS4class{MSMSpepXML}"} is used as
      parameter and a list of matrices or integers of length equal the
      number of spectrum queries is returned.
    }

  }
}

\keyword{methods}
