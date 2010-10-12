\name{spectrumQueries-methods}
\docType{methods}
\alias{spectrumQueries-methods}
\alias{spectrumQueries<-}
\alias{spectrumQueries}

\title{ Returns the spectrum query(ies). }

\description{
  This method returns and sets the (list of) spectrum queries for a
  peptide identificaiton result. The value is set automatically when
  creating an instance of \code{"\linkS4class{MSMSpepXML}"} with
  \code{\link{parseMSMSpepXML}}, and should generally not need to be 
  changed. Note that the \code{[} operator can also be used to extract a
  unique or a list of search hits.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSMSpepXML")}}{
      object of class \code{"\linkS4class{MSMSpepXML}"} is given as
      parameter and a list of \code{"\linkS4class{SpectrumQuery}"}
      instances is returned.
    }
  }
}

\keyword{methods}

