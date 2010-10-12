\name{searchResults-methods}
\docType{methods}
\alias{searchResults-methods}
\alias{searchResults<-}
\alias{searchResults}

\title{ Returns the search result(s) of a spectrum query. }

\description{
  This method returns and sets the (list of) search results for a
  spectrum query object. The value is set automatically when creating an
  instance of \code{"\linkS4class{MSMSpepXML}"} with
  \code{\link{parseMSMSpepXML}}, and should generally not need to be
  changed. Note that the \code{[} operator can also be used to extract a
  unique or a list of search results.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SpectrumQuery")}}{
      object of class \code{"\linkS4class{SpectrumQuery}"} is given as
      parameter and a list of \code{"\linkS4class{SearchResult}"} instances
      is returned.
    }
  }
}
\keyword{methods}

