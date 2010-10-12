\name{proteins-methods}
\docType{methods}
\alias{proteins-methods}
\alias{proteins<-}
\alias{proteins}

\title{ Returns the proteins associated with identifies peptides. }

\description{
  This method returns and sets the protein identifiers for individual
  search hits. The values are set automatically when creating an
  instance of \code{"\linkS4class{MSMSpepXML}"} with
  \code{\link{parseMSMSpepXML}}, and should generally not need to be
  changed.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchHit")}}{
      an object of class \code{"\linkS4class{SearchHit}"} is used as
      parameter and a character vector of possibly multiple proteins is 
      returned. 
    }
    \item{\code{signature(object = "SearchResult")}}{
      an object of class \code{"\linkS4class{SearchResult}"} is used as
      parameter and a character vector of all possibly alternative
      proteins for the search hits is returned.
    }
    \item{\code{signature(object = "SpectrumQuery")}}{
      an object of class \code{"\linkS4class{SpectrumQuery}"} is used as 
      parameter and a character vector of all possibly alternative
      proteins for the search results is returned.
    }
    \item{\code{signature(object = "MSMSpepXML")}}{
      an object of class \code{"\linkS4class{MSMSpepXML}"} is used as 
      parameter and a list of characters of length equal to the number
      of spectrum queries is returned.
    }
  }
}

\keyword{methods}
