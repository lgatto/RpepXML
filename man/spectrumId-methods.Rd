\name{spectrumId-methods}
\docType{methods}
\alias{spectrumId-methods}
\alias{spectrumId<-}
\alias{spectrumId}

\title{ Returns the id of a spectrum query. }

\description{
  This method returns and sets the id of a spectrum query. The value is
  set automatically when creating an instance of
  \code{"\linkS4class{MSMSpepXML}"} with 
  \code{\link{parseMSMSpepXML}}, and should generally not need to be
  changed. 
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SpectrumQuery")}}{
      object of class \code{"\linkS4class{MSMSpepXML}"} is given as
      parameter and an integer is returned.
    }
  }
}
\keyword{methods}

