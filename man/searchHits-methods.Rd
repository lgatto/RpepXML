\name{searchHits-methods}
\docType{methods}
\alias{searchHits-methods}
\alias{searchHits<-}
\alias{searchHits}

\title{ Returns the search hit(s) of a search result. }

\description{
  This method returns and sets the (list of) search hits for a search
  result object. The value is set automatically when creating an
  instance of \code{"\linkS4class{MSMSpepXML}"} with
  \code{\link{parseMSMSpepXML}}, and should generally not need to be
  changed. Note that the \code{[} operator can also be used to extract a
  unique or a list of search hits.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchResult")}}{
      object of class \code{"\linkS4class{SearchResult}"} is given as
      parameter and a list of \code{"\linkS4class{SearchHit}"} instances
      is returned.
    }
  }
}

\keyword{methods}

