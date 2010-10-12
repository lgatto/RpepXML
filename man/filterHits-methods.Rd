\name{filterHits-methods}
\docType{methods}
\alias{filterHits}
\alias{filterHits-methods}

\title{ Method to filter hits based on their rand and ion score.}

\description{
  These methods take two arguments to filter out search hits. When
  applied to an search hit, the method returns a logical indication
  whether the hit passes the filtering criteria. When applied to higher
  level objects, the method returns an object of the same class, where
  all \code{SearchHits} objects that did not pass the filtering have
  been removed. 
}

\usage{ filterHits(object, rank, ionscore)}

\arguments{
  \item{object}{ see object types below. }
  \item{rank}{ integer defining the highest hit rank to be kept,
    default is 1. }
  \item{ionscore}{ numeric definig the lowest ion score to be kept,
    default is 20. }
}


\section{Methods}{ 
  \describe{
    \item{\code{signature(object = "MSMSpepXML")}}{
      an object of class \code{"\linkS4class{MSMSpepXML}"}
    }
    \item{\code{signature(object = "SearchHit")}}{
     an object of class \code{"\linkS4class{SearchHit}"}
    }
    \item{\code{signature(object = "SearchResult")}}{
     an object of class \code{"\linkS4class{SearchResult}"}
    }
    \item{\code{signature(object = "SpectrumQuery")}}{
     an object of class \code{"\linkS4class{SpectrumQuery}"}
    }
  }
}

\keyword{methods}
