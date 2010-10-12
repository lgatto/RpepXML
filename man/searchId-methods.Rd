\name{searchId-methods}
\docType{methods}
\alias{searchId-methods}
\alias{searchId<-}
\alias{searchId}

\title{ Returns the id of the search result. }

\description{
  This method returns and sets the search result. The value is set
  automatically when creating an instance of
  \code{"\linkS4class{MSMSpepXML}"} with 
  \code{\link{parseMSMSpepXML}}, and should generally not need to be
  changed. 
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchResult")}}{
      object of class \code{"\linkS4class{SearchResult}"} is given as
      parameter and an integer is returned.
    }
  }
}
\keyword{methods}

