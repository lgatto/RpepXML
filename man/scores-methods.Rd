\name{scores-methods}
\docType{methods}
\alias{scores-methods}
\alias{scores<-}
\alias{scores}

\title{ Returns the peptide search scores. }

\description{
  This method returns and sets the search scores of a search hit. The
  scores will vary depending on the search engine used. The value is set
  automatically when creating an instance of 
  \code{"\linkS4class{MSMSpepXML}"} with \code{\link{parseMSMSpepXML}},
  and should generally not need to be changed. Mascot's ion score can
  individually be retrieved with \code{\link{ionScore-methods}}. When a
  search result is used, the scores of all the hits are returned.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchHit")}}{
      object of class \code{"\linkS4class{SearchHit}"} is given as
      parameter and a named numeric of length equal to the number of
      scores is returned.
    }
    \item{\code{signature(object = "SearchResult")}}{
      object of class \code{"\linkS4class{SearchResult}"} is given as
      paramter and a matrix is returned where the scores along the rows
      and the search hits along the columns.
    }
  }
}
\keyword{methods}

