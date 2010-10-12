\name{pepFile-methods}
\docType{methods}
\alias{pepFile-methods}
\alias{pepFile}
\alias{pepFile<-}

\title{ Return returns the pepXML file name. }

\description{
  This method returns and sets the full name of the original pepXML file
  used to generate the \code{"\linkS4class{MSMSpepXML}"} object. The value
  is set automatically when creating an instance of
  \code{"\linkS4class{MSMSpepXML}"} with \code{\link{parseMSMSpepXML}},
  and should generally not need to be changed.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSMSpepXML")}}{
      an object of class \code{"\linkS4class{MSMSpepXML}"} is used as
      parameter and a character object is returned.
    }
  }
}
\keyword{methods}

