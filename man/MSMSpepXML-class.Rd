\name{MSMSpepXML-class}
\Rdversion{1.1}
\docType{class}
\alias{MSMSpepXML-class}
\alias{filterHits,MSMSpepXML-method}
\alias{nProteins,MSMSpepXML-method}
\alias{nSearchHits,MSMSpepXML-method}
\alias{nSearchResults,MSMSpepXML-method}
\alias{nSpectrumQueries,MSMSpepXML-method}
\alias{pepFile<-,MSMSpepXML,character-method}
\alias{pepFile<-,MSMSpepXML,ANY-method}
\alias{pepFile,MSMSpepXML-method}
\alias{precNeutralMass,MSMSpepXML-method}
\alias{proteins,MSMSpepXML-method}
\alias{sampleEnzyme<-,MSMSpepXML,character-method}
\alias{sampleEnzyme<-,MSMSpepXML,ANY-method}
\alias{sampleEnzyme,MSMSpepXML-method}
\alias{searchDatabase<-,MSMSpepXML,character-method}
\alias{searchDatabase<-,MSMSpepXML,ANY-method}
\alias{searchDatabase,MSMSpepXML-method}
\alias{searchEngine<-,MSMSpepXML,character-method}
\alias{searchEngine<-,MSMSpepXML,ANY-method}
\alias{searchEngine,MSMSpepXML-method}
\alias{show,MSMSpepXML-method}
\alias{spectrumQueries<-,MSMSpepXML,list-method}
\alias{spectrumQueries<-,MSMSpepXML,ANY-method}
\alias{[,MSMSpepXML-method}
\alias{spectrumQueries,MSMSpepXML-method}
\alias{assumedCharge,MSMSpepXML-method}
\alias{filterHits,MSMSpepXML-method}
\alias{pepSequence,MSMSpepXML-method}
\alias{scores,MSMSpepXML-method}


\title{
  Class to contain MSMS peptide identification data.
}

\description{Class for the \code{msms_pipeline_analysis} element of the
  pepXML standard. This is the main object that will hold all the MSMS
  database search results.
}

\section{Slots}{
  \describe{
    \item{\code{pepFile}:}{Object of class \code{"character"} containing
  the name of the file used parsed to create the \code{MSMSpepXML} object.}
    \item{\code{searchEngine}:}{Object of class \code{"character"}
  describes the search engine used to generate data.}
    \item{\code{sampleEnzyme}:}{Object of class \code{"character"}
  described the enzyme used for sample digestion.}
    \item{\code{searchDatabase}:}{Object of class \code{"character"}
  describing the data base used for sepctra identification.}
\item{\code{spectrumQueries}:}{Object of class \code{"list"}
      containing identification results for the queried MSMS sepctra. 
      See also \code{"\linkS4class{SpectrumQuery}"}.}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"},
  documenting the version of the class.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  \describe{
    \item{filterHits}{\code{signature(object = "MSMSpepXML")}: see \code{\link{filterHits-methods}}.}
    \item{nProteins}{\code{signature(object = "MSMSpepXML")}: see \code{\link{nProteins-methods}}.}
    \item{nSearchHits}{\code{signature(object = "MSMSpepXML")}: see \code{\link{nSearchHits-methods}}.}
    \item{nSearchResults}{\code{signature(object = "MSMSpepXML")}: see \code{\link{nSearchResults-methods}}.}
    \item{pepFile<-}{\code{signature(object = "MSMSpepXML", value = "character")}: see \code{\link{pepFile-methods}}.}
    \item{pepFile}{\code{signature(object = "MSMSpepXML")}: see \code{\link{pepFile-methods}}.}
    \item{proteins}{\code{signature(object = "MSMSpepXML")}: see \code{\link{proteins-methods}}.}
    \item{sampleEnzyme<-}{\code{signature(object = "MSMSpepXML", value = "character")}: see \code{\link{sampleEnzyme-methods}}.}
    \item{sampleEnzyme}{\code{signature(object = "MSMSpepXML")}: see \code{\link{sampleEnzyme-methods}}.}
    \item{searchDatabase<-}{\code{signature(object = "MSMSpepXML", value = "character")}: see \code{\link{searchDatabase-methods}}.}
    \item{searchDatabase}{\code{signature(object = "MSMSpepXML")}: see \code{\link{searchDatabase-methods}}.}
    \item{searchEngine<-}{\code{signature(object = "MSMSpepXML", value = "character")}: see \code{\link{searchEngine-methods}}.}
    \item{searchEngine}{\code{signature(object = "MSMSpepXML")}: see \code{\link{searchEngine-methods}}.}
    \item{show}{\code{signature(object = "MSMSpepXML")}: see \code{\link{show-methods}}.}
    \item{spectrumQueries<-}{\code{signature(object = "MSMSpepXML", value = "list")}: see \code{\link{spectrumQueries-methods}}.}
    \item{spectrumQueries}{\code{signature(object = "MSMSpepXML")}: see \code{\link{spectrumQueries-methods}}.}
	 }
}

\author{
  Laurent Gatto
}

\seealso{
  \code{\link{parseMSMSpepXML}},
  \code{"\linkS4class{SpectrumQuery}"},
  \code{"\linkS4class{SearchResult}"},
  \code{"\linkS4class{SearchHit}"}
}
\examples{
showClass("MSMSpepXML")
}
\keyword{classes}
