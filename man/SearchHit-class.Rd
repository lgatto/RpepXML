\name{SearchHit-class}
\Rdversion{1.1}
\docType{class}
\alias{SearchHit-class}
\alias{filterHits,SearchHit-method}
\alias{hitRank<-,SearchHit,integer-method}
\alias{hitRank<-,SearchHit,ANY-method}
\alias{hitRank,SearchHit-method}
\alias{ionScore,SearchHit-method}
\alias{modifications<-,SearchHit,data.frame-method}
\alias{modifications<-,SearchHit,ANY-method}
\alias{modifications,SearchHit-method}
\alias{nProteins,SearchHit-method}
\alias{pepSequence<-,SearchHit,character-method}
\alias{pepSequence<-,SearchHit,ANY-method}
\alias{pepSequence,SearchHit-method}
\alias{proteins<-,SearchHit,character-method}
\alias{proteins<-,SearchHit,ANY-method}
\alias{proteins,SearchHit-method}
\alias{scores<-,SearchHit,numeric-method}
\alias{scores<-,SearchHit,ANY-method}
\alias{scores,SearchHit-method}
\alias{show,SearchHit-method}

\title{
  Class to contain individual search hit identification data.
}

\description{
  Class implementing the \code{search_hit} element of the pepXML standard. 
}

\section{Slots}{
  \describe{
    \item{\code{hitRank}:}{Object of class \code{"integer"}, indicating
      the hit rank of the current search hit. }
    \item{\code{proteins}:}{Object of class \code{"character"},
      documenting the protein(s) in which the current peptide hit is
      found. }
    \item{\code{pepSequence}:}{Object of class \code{"character"}, the
      sequence of the peptide. }
    \item{\code{modifications}:}{Object of class \code{"data.frame"},
      documenting the modifications identified in the current peptide. }
    \item{\code{scores}:}{Object of class (names) \code{"numeric"}, documenting
      the different search scores values. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"},
      documenting the class version. }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  \describe{
    \item{filterHits}{\code{signature(object = "SearchHit")}: see \code{\link{filterHits-methods}}. }
    \item{hitRank<-}{\code{signature(object = "SearchHit", value = "integer")}: see \code{\link{hitRank-methods}}. }
    \item{hitRank}{\code{signature(object = "SearchHit")}: see \code{\link{hitRank-methods}}. }
    \item{ionScore}{\code{signature(object = "SearchHit")}: see \code{\link{ionScore-methods}}. }
    \item{modifications<-}{\code{signature(object = "SearchHit", value = "data.frame")}: see \code{\link{modifications-methods}}. }
    \item{modifications}{\code{signature(object = "SearchHit")}: see \code{\link{modifications-methods}}. }
    \item{nProteins}{\code{signature(object = "SearchHit")}: see \code{\link{nProteins-methods}}. }
    \item{pepSequence<-}{\code{signature(object = "SearchHit", value = "character")}: see \code{\link{pepSequence-methods}}. }
    \item{pepSequence}{\code{signature(object = "SearchHit")}: see \code{\link{pepSequence-methods}}. }
    \item{proteins<-}{\code{signature(object = "SearchHit", value = "character")}: see \code{\link{proteins-methods}}. }
    \item{proteins}{\code{signature(object = "SearchHit")}: see \code{\link{proteins-methods}}. }
    \item{scores<-}{\code{signature(object = "SearchHit", value = "numeric")}: see \code{\link{scores-methods}}. }
    \item{scores}{\code{signature(object = "SearchHit")}: see \code{\link{scores-methods}}. }
    \item{show}{\code{signature(object = "SearchHit")}: see \code{\link{show-methods}}. }
	 }
}

\author{
  Laurent Gatt
}

\seealso{
  \code{"\linkS4class{MSMSpepXML}"},
  \code{"\linkS4class{SpectrumQuery}"},
  \code{"\linkS4class{SearchResult}"}
}

\examples{
showClass("SearchHit")
}
\keyword{classes}
