\name{proteinSummary}
\alias{proteinSummary}

\title{
  Return a table of identified proteins.
}

\description{
  This function works with instances of
  \code{"\linkS4class{MSMSpepXML}"}, 
  \code{"\linkS4class{SpectrumQuery}"},
  \code{"\linkS4class{SearchResult}"} or
  \code{"\linkS4class{SearchHit}"} classes and returns a table summary
  of the proteins inferred from the peptide hits. All proteins,
  irrespective of peptide scores of hit rank are return.
}

\usage{
proteinSummary(object)
}


\arguments{
  \item{object}{
    object of class \code{"\linkS4class{MSMSpepXML}"}, 
  \code{"\linkS4class{SpectrumQuery}"},
  \code{"\linkS4class{SearchResult}"} or
  \code{"\linkS4class{SearchHit}"}.
  }
}

\value{
  A table.
}

\author{
  Laurent Gatto
}

\seealso{
  See Also \code{\link{proteins-methods}}.
}

\examples{
 pepxml <- system.file("extdata","tinySearch.pepXML",package="RpepXML")
 x <- parseMSMSpepXML(pepxml)
 proteinSummary(x)
}

\keyword{ utilities }
\keyword{ manip }
