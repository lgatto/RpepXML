\name{parseMSMSpepXML}
\alias{parseMSMSpepXML}

\title{
  Read a pepXML file into R
}

\description{
  This function reads a pepXML given as parameter and parses it,
  producing an object of class \code{"\linkS4class{MSMSpepXML}"}.
}

\usage{
parseMSMSpepXML(pepfile, verbose = TRUE)
}


\arguments{
  \item{pepfile}{
    the name of the pepXML file which the data are to be read from.
  }
  \item{verbose}{
    a logical value indicating if a progress bar should be displayed.
  }
}

\value{
  an object of class \code{"\linkS4class{MSMSpepXML}"}.
}


\author{
Laurent Gatto
}

\note{
  The pepXML file should not contain any XML comments as these are not
  (yet) handled properly and will be considered as nodes in the XML
  tree. 
}


\examples{
 pepxml <- system.file("extdata","tinySearch.pepXML",package="RpepXML")
 x <- parseMSMSpepXML(pepxml)
 show(x)
}


\keyword{ utilities }
\keyword{ manip }

