\name{searchEngine-methods}
\docType{methods}
\alias{searchEngine-methods}
\alias{searchEngine<-}
\alias{searchEngine}

\title{ Returns the name of the search engine. }

\description{
  This method returns and sets the search engine name as recorded in the 
  original pepXML file. The value is set automatically when creating an
  instance of \code{"\linkS4class{MSMSpepXML}"} with
  \code{\link{parseMSMSpepXML}}, and should generally not need to be
  changed. 
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSMSpepXML")}}{
      object of class \code{"\linkS4class{MSMSpepXML}"} is given as
      parameter and a character is returned.
    }
  }
}
\keyword{methods}

