\name{hitRank-methods}
\docType{methods}
\alias{hitRank-methods}
\alias{hitRank<-}
\alias{hitRank}

\title{ Accessor methods of the \code{hitRank} slot }

\description{
  These methods allow to retrieve the \code{hitRank} slot
  of instances of all the instance of \code{"\linkS4class{SearchHit}"}
  from an \code{"\linkS4class{SearchResult}"} or individual
  \code{"\linkS4class{SearchHit}"} and update the slot of the
  latter. The values are set automatically when creating an instance of 
  \code{"\linkS4class{MSMSpepXML}"} with \code{\link{parseMSMSpepXML}},
  and should generally not need to be changed.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchHit")}}{
      object of class \code{"\linkS4class{SearchHit}"}
    }
    \item{\code{signature(object = "SearchResult")}}{
      object of class \code{"\linkS4class{SearchResult}"}
    }
  }
}

\keyword{methods}
