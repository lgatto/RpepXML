\name{SpectrumQuery-class}
\Rdversion{1.1}
\docType{class}
\alias{SpectrumQuery-class}
\alias{assumedCharge<-,SpectrumQuery,numeric-method}
\alias{assumedCharge<-,SpectrumQuery,ANY-method}
\alias{assumedCharge,SpectrumQuery-method}
\alias{filterHits,SpectrumQuery-method}
\alias{nProteins,SpectrumQuery-method}
\alias{nSearchHits,SpectrumQuery-method}
\alias{nSearchResults,SpectrumQuery-method}
\alias{precNeutralMass<-,SpectrumQuery,numeric-method}
\alias{precNeutralMass<-,SpectrumQuery,ANY-method}
\alias{precNeutralMass,SpectrumQuery-method}
\alias{proteins,SpectrumQuery-method}
\alias{queryIndex<-,SpectrumQuery,integer-method}
\alias{queryIndex<-,SpectrumQuery,ANY-method}
\alias{queryIndex,SpectrumQuery-method}
\alias{searchResults<-,SpectrumQuery,list-method}
\alias{searchResults<-,SpectrumQuery,ANY-method}
\alias{[,SpectrumQuery-method}
\alias{searchResults,SpectrumQuery-method}
\alias{show,SpectrumQuery-method}
\alias{spectrumId<-,SpectrumQuery,character-method}
\alias{spectrumId<-,SpectrumQuery,ANY-method}
\alias{spectrumId,SpectrumQuery-method}
\alias{startScan,SpectrumQuery-method}
\alias{startScan<-,SpectrumQuery,integer-method}
\alias{startScan<-,SpectrumQuery,ANY-method}
\alias{endScan,SpectrumQuery-method}
\alias{endScan<-,SpectrumQuery,integer-method}
\alias{endScan<-,SpectrumQuery,ANY-method}

\title{
  Class to contain individual spectrum query identification data.
}

\description{
  Class implementing the \code{spectrum_query} element of the pepXML standard. 
}

\section{Slots}{
  \describe{
    \item{\code{searchResults}:}{Object of class \code{"list"}
  (generally of length 1), containing the search results, see \code{"\linkS4class{SearchResult}"}. }
    \item{\code{spectrumId}:}{Object of class \code{"character"}, id of
  the spectrum. }
    \item{\code{startScan}:}{Object of class \code{"integer"} providing
  the start scan. }
    \item{\code{endScan}:}{Object of class \code{"integer"} providing
  the end scan. }
    \item{\code{precNeutralMass}:}{Object of class \code{"numeric"}
  giving the precursor's neutral mass. }
    \item{\code{assumedCharge}:}{Object of class \code{"numeric"}
  documenting the precursor's assumed charge. }
    \item{\code{queryIndex}:}{Object of class \code{"integer"}
  specifying the query unique index. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"},
  documenting the version of the class.}
  }
}
\section{Extends}{
Class \code{\link[Biobase]{Versioned}}, directly.
}
\section{Methods}{
  \describe{
    \item{assumedCharge<-}{\code{signature(object = "SpectrumQuery", value = "numeric")}: ... }
    \item{assumedCharge}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{endScan<-}{\code{signature(object = "SpectrumQuery", value = "integer")}: ... }
    \item{endScan}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{filterHits}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{nProteins}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{nSearchHits}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{nSearchResults}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{precNeutralMass<-}{\code{signature(object = "SpectrumQuery", value = "numeric")}: ... }
    \item{precNeutralMass}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{proteins}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{queryIndex<-}{\code{signature(object = "SpectrumQuery", value = "integer")}: ... }
    \item{queryIndex}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{searchResults<-}{\code{signature(object = "SpectrumQuery", value = "list")}: ... }
    \item{searchResults}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{show}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{spectrumId<-}{\code{signature(object = "SpectrumQuery", value = "character")}: ... }
    \item{spectrumId}{\code{signature(object = "SpectrumQuery")}: ... }
    \item{startScan<-}{\code{signature(object = "SpectrumQuery", value = "integer")}: ... }
    \item{startScan}{\code{signature(object = "SpectrumQuery")}: ... }
	 }
}

\author{
  Laurent Gatto
}

\seealso{
  \code{"\linkS4class{MSMSpepXML}"},
  \code{"\linkS4class{SearchResult}"},
  \code{"\linkS4class{SearchHit}"}
}
\examples{
showClass("SpectrumQuery")
}
\keyword{classes}
