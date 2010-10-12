\name{SearchResult-class}
\Rdversion{1.1}
\docType{class}
\alias{SearchResult-class}
\alias{filterHits,SearchResult-method}
\alias{hitRank,SearchResult-method}
\alias{nProteins,SearchResult-method}
\alias{nSearchHits,SearchResult-method}
\alias{pepSequence,SearchResult-method}
\alias{proteins,SearchResult-method}
\alias{scores,SearchResult-method}
\alias{searchHits<-,SearchResult,list-method}
\alias{searchHits<-,SearchResult,ANY-method}
\alias{searchHits,SearchResult-method}
\alias{[,SearchResult-method}
\alias{searchId<-,SearchResult,integer-method}
\alias{searchId<-,SearchResult,ANY-method}
\alias{searchId,SearchResult-method}
\alias{show,SearchResult-method}

\title{
  Class to contain individual search result data.
}

\description{
  Class implementing the \code{search_result} of the pepXML standard. 
}

\section{Slots}{
  \describe{
    \item{\code{searchHits}:}{Object of class \code{"list"} containing a
      list of search hits, see \code{"\linkS4class{SearchHit}"}.}
    \item{\code{searchId}:}{Object of class \code{"integer"}, id of the
      search hit. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"},
      documenting the version of the class. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  \describe{
    \item{filterHits}{\code{signature(object = "SearchResult")}: see \code{\link{filterHits-methods}}. }
    \item{hitRank}{\code{signature(object = "SearchResult")}: see \code{\link{hitRank-methods}}. }
    \item{nProteins}{\code{signature(object = "SearchResult")}: see \code{\link{nProteins-methods}}. }
    \item{nSearchHits}{\code{signature(object = "SearchResult")}: see \code{\link{nSearchHits-methods}}. }
    \item{pepSequence}{\code{signature(object = "SearchResult")}: see \code{\link{pepSequence-methods}}. }
    \item{proteins}{\code{signature(object = "SearchResult")}: see \code{\link{proteins-methods}}. }
    \item{scores}{\code{signature(object = "SearchResult")}: see \code{\link{scores-methods}}. }
    \item{searchHits<-}{\code{signature(object = "SearchResult", value = "list")}: see \code{\link{searchHits-methods}}. }
    \item{searchHits}{\code{signature(object = "SearchResult")}: see \code{\link{searchHits-methods}}. }
    \item{searchId<-}{\code{signature(object = "SearchResult", value = "integer")}: see \code{\link{searchId-methods}}. }
    \item{searchId}{\code{signature(object = "SearchResult")}: see \code{\link{searchId-methods}}. }
    \item{show}{\code{signature(object = "SearchResult")}: see \code{\link{show-methods}}. }
  }
}

\author{
  Laurent Gatto
}

\seealso{
  \code{"\linkS4class{MSMSpepXML}"},
  \code{"\linkS4class{SpectrumQuery}"},
  \code{"\linkS4class{SearchHit}"}
}
\examples{
showClass("SearchResult")
}
\keyword{classes}
