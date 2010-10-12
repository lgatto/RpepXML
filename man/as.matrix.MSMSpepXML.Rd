\name{as.matrix.MSMSpepXML}
\alias{as.matrix.MSMSpepXML}

\title{
  Transforms an \code{MSMSpepXML} object in a simplified matrix.
}

\description{
  This function transforms an \code{"\linkS4class{MSMSpepXML}"} instance
  in a matrix. The object is significantly simplified, as only the first
  search results and hits are kept. The matrix includes query index,
  spectrum id, ion score, peptide sequence, the number of proteins 
  inferred and the protein(s) name(s).
}

\usage{
as.matrix.MSMSpepXML(x, ..., all = TRUE)
}

\arguments{
  \item{x}{
    an object of class \code{"\linkS4class{MSMSpepXML}"}.
  }
  \item{...}{
    further options passed to next methods.
  }
  \item{all}{
    logical, whether all rows or only the ones with identified peptides
    should ne returned.
  }
}

\value{
  A matrix.
}

\author{
  Laurent Gatto
}


\section{Warning }{
  If search results and hits are removed, warnings are issued.
}

\examples{
 pepxml <- system.file("extdata","tinySearch.pepXML",package="RpepXML")
 x <- parseMSMSpepXML(pepxml)
 as.matrix(x)
 as.matrix(filterHits(x),all=TRUE)
 as.matrix(filterHits(x),all=FALSE)
}

\keyword{ utilities }
\keyword{ manip }
