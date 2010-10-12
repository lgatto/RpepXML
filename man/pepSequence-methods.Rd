\name{pepSequence-methods}
\docType{methods}
\alias{pepSequence-methods}
\alias{pepSequence<-}
\alias{pepSequence}

\title{ Returns the sequence(s) of an identified peptide. }

\description{
  This function returns and sets the sequence of the peptide(s)
  identified for one search hit or all the search hits of a search
  result. The output is a character of length equal to the number of
  search hits (when quering a search result) where each returned string
  has the format 'previous amino accid'.'peptide sequence'.'next amino
  acid'. A '-' represents the end of the protein sequence. The values
  are set automatically when creating an instance of
  \code{"\linkS4class{MSMSpepXML}"} with \code{\link{parseMSMSpepXML}}, 
  and should generally not need to be changed.

}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "SearchHit")}}{
      an object of class \code{"\linkS4class{SearchHit}"} is used as
      parameter and a character is returned.      
    }
    \item{\code{signature(object = "SearchResult")}}{
      an object of class \code{"\linkS4class{SearchResult}"} is used as
      parameter and a character vector of length equal to the number of
      search hits is returned.
    }
  }
}

\keyword{methods}
