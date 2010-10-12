\name{RpepXML-package}
\alias{RpepXML-package}
\alias{RpepXML}
\docType{package}
\title{R interface to the pepXML standard}
\description{
  RpepXML allows to load MSMS pipeline analysis 
  data in pepXML format and provides programmatic
  access to the different elements. It also allow to
  filter the identificaiton results based on peptide rank
  and ion score.
}
\details{
  \tabular{ll}{
    Package: \tab RpepXML\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2010-08-27\cr
    License: \tab Artitic-2.0\cr
    LazyLoad: \tab yes\cr
  }

  From the Trans-Proteomic Pipeline (TPP) pepXML page:
  \sQuote{pepXML Is an open data format developed at the
    SPC/Institute for Systems biology for the storage,
    exchange, and processing of peptide sequence assignments
    of MS/MS scans. pepXML is intended to provide a common
    data output format for many different MS/MS search
    engines and subsequent peptide-level analyses.}

  Currently, only the \code{msms_pipeline_analysis} element and children
  are implemented. The overall structure is simplified and not all
  attributes are implemented (yet). The most important elements should
  be in place to import, manipulate and filter pepXML data.
  
  See the package vignette for details about the package's usage
  and underlying classes. An dummy example data set
  \code{tinySearch.pepXML} is also provided.
}

\author{
Laurent Gatto
Maintainer: Laurent Gatto <lg390@cam.ac.uk>
}

\references{
  \url{http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML},
  \url{http://sashimi.sourceforge.net/schema_revision/pepXML/Docs/pepXML_v18.html}.
}

\keyword{ package }
