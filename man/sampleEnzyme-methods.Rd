\name{sampleEnzyme-methods}
\docType{methods}
\alias{sampleEnzyme-methods}
\alias{sampleEnzyme<-}
\alias{sampleEnzyme}


\title{ Returns the sample enzyme. }

\description{
  This method returns and sets the sample enzyme name as recorded in the 
  original pepXML file The value is set automatically when creating an
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

