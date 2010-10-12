\name{queryIndex-methods}
\docType{methods}
\alias{queryIndex-methods}
\alias{queryIndex<-}
\alias{queryIndex}

\title{ Returns the spectrum query index. }
\description{
  This method retrieves and sets the spectrum index number of a spectrum
  query. This value is set automatically when creating an instance of
  \code{"\linkS4class{MSMSpepXML}"} with \code{\link{parseMSMSpepXML}},
  and should generally not need to be changed.
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "SpectrumQuery")}}{
      object of class \code{"\linkS4class{SpectrumQuery}"} is given as
      parameter and an integer is returned.}
  }
}

\keyword{methods}

