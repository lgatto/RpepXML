parseMSMSpepXML <- function(pepfile,verbose=TRUE) {
  ## =========================================
  ## Parameters:
  ##   file: a pep.xml file
  ##   verbose: logical, whether progress bar is shown
  ## Return:
  ##   an object of class MSMSpepXML
  ## -----------------------------------------
  ifelse(verbose,progress <- "text",progress <- "none")
  doc <- xmlRoot(xmlTreeParse(pepfile))
  msms_run_summary <- doc[[1]]
  ## New instance of MSMSpepXML
  ret <- new("MSMSpepXML")
  ## Populate MSMSpepXML's attributes
  pepFile(ret) <- pepfile
  sampleEnzyme(ret) <- xmlGetAttr(msms_run_summary[["sample_enzyme"]],"name")
  searchEngine(ret) <- xmlGetAttr(msms_run_summary[["search_summary"]],"search_engine")
  searchDatabase(ret) <- xmlGetAttr(msms_run_summary[["search_summary"]][["search_database"]],"database_release_identifier")
  n <- which(names(msms_run_summary)=="spectrum_query")
  spqrieslist <- msms_run_summary[n]
  ## Generating a list of SpectrumQuery objects
  spectrumQueries(ret) <- llply(spqrieslist,function(sq) {
    ## New instance of SpectrumQuery object
    x <- new("SpectrumQuery")
    ## Populate SpectrumQuery's attributes
    spectrumId(x) <- xmlGetAttr(sq,"spectrum")
    startScan(x) <- as.integer(xmlGetAttr(sq,"start_scan"))
    endScan(x) <- as.integer(xmlGetAttr(sq,"end_scan"))
    precNeutralMass(x) <- as.numeric(xmlGetAttr(sq,"precursor_neutral_mass"))
    assumedCharge(x) <- as.numeric(xmlGetAttr(sq,"assumed_charge"))
    queryIndex(x) <- as.integer(xmlGetAttr(sq,"index"))
    ## List of SearchResult objects
    searchResults(x) <- list()
    if (length(sq)>0)
      searchResults(x) <- getSearchResults(sq)
    return(x)
  },.progress=progress)
  return(ret)
}

getSearchResults <- function(spqry) {
  ## =========================================
  ## Parameters:
  ##   spqry: xml entry for a spectrum_query
  ## Return:
  ##   a list of SearchResult objects
  ## -----------------------------------------
  if (length(spqry)>0) {
    return(xmlApply(spqry,function(x) {
      y <- new("SearchResult")
      searchId(y) <- as.integer(xmlGetAttr(x,"search_id"))
      searchHits(y) <- getSearchHits(x)
      return(y)
    }))
  } else {
    return(list)
    }
}

getSearchHits <- function(srchres) {
  ## =========================================
  ## Parameters:
  ##   srchres: xml entry for a search_result
  ## Return:
  ##   a list of SearchHit objects
  ## -----------------------------------------
  return(xmlApply(srchres, function(x) {
    y <- new("SearchHit")
    hitRank(y) <- as.integer(xmlGetAttr(x,"hit_rank"))
    pepSequence(y) <- paste(xmlGetAttr(x,"peptide_prev_aa"),
                            xmlGetAttr(x,"peptide"),
                            xmlGetAttr(x,"peptide_next_aa"),
                            sep=".")
    proteins(y) <- c(xmlGetAttr(x,"protein"),
                     unlist(xmlSApply(x,xmlGetAttr,"protein")))
    scores(y) <- as.numeric(unlist(xmlSApply(x,xmlGetAttr,"value")))
    names(scores(y)) <- unlist(xmlSApply(x,xmlGetAttr,"name"))
    if (any(names(x)=="modification_info")) {
      modifications(y) <- as.data.frame(xmlSApply(x[["modification_info"]],xmlAttrs))
    }
    return(y)
  }))
}


parseMSMSpepXMLfile.old <- function(pepfile,verbose=TRUE) {
  ## =========================================
  ## Parameters:
  ##   file: a pep.xml file
  ##   verbose: logical, whether progress bar is shown
  ## Return:
  ##   an object of class MSMSpepXML
  ## -----------------------------------------  
  doc <- xmlRoot(xmlTreeParse(pepfile))
  msms_run_summary <- doc[[1]]
  ## New instance of MSMSpepXML
  ret <- new("MSMSpepXML")
  ## Populate MSMSpepXML's attributes
  pepFile(ret) <- pepfile
  sampleEnzyme(ret) <- xmlGetAttr(msms_run_summary[["sample_enzyme"]],"name")
  searchEngine(ret) <- xmlGetAttr(msms_run_summary[["search_summary"]],"search_engine")
  searchDatabase(ret) <- xmlGetAttr(msms_run_summary[["search_summary"]][["search_database"]],"database_release_identifier")
  n <- which(names(msms_run_summary)=="spectrum_query")
  ## List of SpectrumQuery objects
  spqries <- vector("list",length=length(n))
  if (verbose)
    pb <- txtProgressBar(min=0,max=length(n),style=3)
  ## Populating MSMSpepXML's spectrum queries
  for (i in n) {
    ## New instance of SpectrumQuery object
    x <- new("SpectrumQuery")
    ## Populate SpectrumQuery's attributes
    spectrumId(x) <- xmlGetAttr(msms_run_summary[[i]],"spectrum")
    startScan(x) <- as.integer(xmlGetAttr(msms_run_summary[[i]],"start_scan"))
    endScan(x) <- as.integer(xmlGetAttr(msms_run_summary[[i]],"end_scan"))
    precNeutralMass(x) <- as.numeric(xmlGetAttr(msms_run_summary[[i]],"precursor_neutral_mass"))
    assumedCharge(x) <- as.numeric(xmlGetAttr(msms_run_summary[[i]],"assumed_charge"))
    queryIndex(x) <- as.integer(xmlGetAttr(msms_run_summary[[i]],"index"))
    ## List of SearchResult objects
    searchResults(x) <- list()
    if (length(msms_run_summary[[i]])>0)
      searchResults(x) <- getSearchResults.old(msms_run_summary[[i]])
    spqries[[i-2]] <- x
    if (verbose)
      setTxtProgressBar(pb, i)
  }
  if (verbose)
    close(pb)
  spectrumQueries(ret) <- spqries
  return(ret)
}

getSearchResults.old <- function(spqry) {
  ## =========================================
  ## Parameters:
  ##   spqry: xml entry for a spectrum_query
  ## Return:
  ##   a list of SearchResult objects
  ## -----------------------------------------
  ret <- vector(mode="list",length=length(spqry))
  for (i in 1:length(spqry)) {
    x <- new("SearchResult")
    searchId(x) <- as.integer(xmlGetAttr(spqry[[i]],"search_id"))
    searchHits(x) <- getSearchHits.old(spqry[[i]])
    ret[[i]] <- x
  }
  return(ret)
}

getSearchHits.old <- function(srchres) {
  ## =========================================
  ## Parameters:
  ##   srchres: xml entry for a search_result
  ## Return:
  ##   a list of SearchHit objects
  ## -----------------------------------------
  ret <- vector(mode="list",length=length(srchres))
  for (i in 1:length(srchres)) {
    x <- new("SearchHit")
    hitRank(x) <- as.integer(xmlGetAttr(srchres[[i]],"hit_rank"))
    pepSequence(x) <- paste(xmlGetAttr(srchres[[i]],"peptide_prev_aa"),
                            xmlGetAttr(srchres[[i]],"peptide"),
                            xmlGetAttr(srchres[[i]],"peptide_next_aa"),
                            sep=".")
    proteins(x) <- c(xmlGetAttr(srchres[[i]],"protein"),
                     unlist(xmlSApply(srchres[[i]],xmlGetAttr,"protein")))
    scores(x) <- as.numeric(unlist(xmlSApply(srchres[[i]],xmlGetAttr,"value")))
    names(scores(x)) <- unlist(xmlSApply(srchres[[i]],xmlGetAttr,"name"))
    if (any(names(srchres[[i]])=="modification_info")) {
      modifications(x) <- as.data.frame(xmlSApply(srchres[[i]][["modification_info"]],xmlAttrs))
    }
    ret[[i]] <- x
  }
  return(ret)
}
