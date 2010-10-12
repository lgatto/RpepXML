\name{startScan-methods}
\docType{methods}
\alias{startScan-methods}
\alias{endScan-methods}
\alias{startScan<-}
\alias{startScan}
\alias{endScan<-}
\alias{endScan}



\title{ Accessor of the \code{startScan} and \code{endScan} slots }

\description{
  These methods allow to retrieve or update the \code{startScan} and
  \code{endScan} slots of an instance of
  \code{"\linkS4class{SearchResult}"}. The values are set automatically
  when creating an instance of \code{"\linkS4class{MSMSpepXML}"} with
  \code{\link{parseMSMSpepXML}}, and should generally not need to be
  changed. 
}

\section{Methods}{
\describe{

  \item{\code{signature(object = "SpectrumQuery")}}{
    object of class \code{"\linkS4class{SpectrumQuery}"}
  }
}}

\keyword{methods}
