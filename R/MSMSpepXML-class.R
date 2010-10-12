### =================================================
## Upper-level class for a simplified element
## 'msms_pipeline_analysis' of the pepXML standard.
## See pepXML specifications for underlying xsd
## -------------------------------------------------

setClass("MSMSpepXML",
         representation(pepFile="character",
                        searchEngine="character",
                        sampleEnzyme="character",
                        searchDatabase="character",
                        spectrumQueries="list"), ## of SpectrumQuery objects
         contains = c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(RpepXML="0.0.1")))
         )        

## ============================
## Accessor method definitions
## ----------------------------
setMethod("pepFile","MSMSpepXML",function(object) object@pepFile)
setMethod("pepFile<-","MSMSpepXML",
          function(object,value="character") object@pepFile <- value)
setReplaceMethod("pepFile",
                 signature(object="MSMSpepXML",
                           value="character"),
                 function(object,value) {
                   object@pepFile <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("searchEngine","MSMSpepXML",function(object) object@searchEngine)
setMethod("searchEngine<-","MSMSpepXML",
          function(object,value="character") object@searchEngine <- value)
setReplaceMethod("searchEngine",
                 signature(object="MSMSpepXML",
                           value="character"),
                 function (object,value) {
                   object@searchEngine <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("sampleEnzyme","MSMSpepXML",function(object) object@sampleEnzyme)
setMethod("sampleEnzyme<-","MSMSpepXML",
          function(object,value="character") object@sampleEnzyme <- value)
setReplaceMethod("sampleEnzyme",
                 signature(object="MSMSpepXML",
                           value="character"),
                 function (object,value) {
                   object@sampleEnzyme <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("searchDatabase","MSMSpepXML",function(object) object@searchDatabase)
setMethod("searchDatabase<-","MSMSpepXML",
          function(object,value="character") object@searchDatabase <- value)
setReplaceMethod("searchDatabase",
                 signature(object="MSMSpepXML",
                           value="character"),
                 function (object,value) {
                   object@searchDatabase <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("spectrumQueries","MSMSpepXML",function(object) object@spectrumQueries)
setMethod("spectrumQueries<-","MSMSpepXML",
          function(object,value="list") object@spectrumQueries <- value)
setReplaceMethod("spectrumQueries",
                 signature(object="MSMSpepXML",
                           value="list"),
                 function (object,value) {
                   object@spectrumQueries <- value
                   if (validObject(object))
                     return(object)
                 })

## =========================
## Other method definitions
## -------------------------
setMethod("show","MSMSpepXML",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat(" File:",basename(pepFile(object)),"\n")
            cat(" Database:",searchDatabase(object),"\n")
            cat(" Engine:",searchEngine(object),"\n")            
            cat(" Enzyme:",sampleEnzyme(object),"\n")
            cat("",nSpectrumQueries(object),"spectrum queries\n")
            cat("",sum(nSearchResults(object)),"peptides identified\n")
          })


setMethod("nSearchResults","MSMSpepXML",          
          function(object) unlist(sapply(spectrumQueries(object),nSearchResults)))

setMethod("nSearchHits","MSMSpepXML",          
          function(object) sapply(spectrumQueries(object),nSearchHits))

setMethod("nSpectrumQueries","MSMSpepXML",
          function(object) return(length(spectrumQueries(object))))

setMethod("proteins","MSMSpepXML",
          function(object) {
            sapply(spectrumQueries(object),proteins)
          })

setMethod("assumedCharge","MSMSpepXML",
          function(object) {
            sapply(spectrumQueries(object),assumedCharge)
          })

setMethod("filterHits","MSMSpepXML",
          function(object,rank=1,ionscore=20) {
            spectrumQueries(object) <-
              lapply(spectrumQueries(object),filterHits,rank,ionscore)
            return(object)
          })



proteinSummary <- function(object) {
  table(unlist(proteins(object)))
}

setMethod("nProteins","MSMSpepXML",
          function(object) lapply(spectrumQueries(object),nProteins))


as.matrix.MSMSpepXML <- function(x,...,all=TRUE) {
  ## get the index of an spectrum query that has
  ## search hits
  i <- which(sapply(nSearchHits(x),function(x) x[1]>0))[1]
  if (is.na(i)) {
    warning("No search hits found.")
    return(NULL)
  }
  scorenames <- names(scores(x[i][1][1]))
  ret <- lapply(spectrumQueries(x),
                function(y) {
                  rw <- c(queryIndex(y),spectrumId(y))
                  if (nSearchResults(y)==0) {
                    rw <- c(rw,"",rep("",length(scorenames)),"0","")
                    return(rw)
                  } else {
                    if (nSearchResults(y)>1) 
                      warning("Only keeping first search result at query index ",queryIndex(y),call.=FALSE)
                    if (nSearchHits(y)[1]>1)
                      warning("Only keeping first search hit at query index ",queryIndex(y),call.=FALSE)
                    sh <- searchHits(searchResults(y)[[1]])[[1]]
                    rw <- c(rw,scores(sh),pepSequence(sh),nProteins(sh),paste(proteins(sh),collapse=";"))
                    return(rw)
                  }
                })
  ret <- t(as.data.frame(ret))
  colnames(ret) <- c("QueryIndex","SpectrumId",scorenames,
                     "Sequence","NumProteins","Proteins")
  rownames(ret) <- 1:nrow(ret)
  if (!all) {
    sel <- ret[,"NumProteins"]!="0"
    ret <- ret[sel,]
  }
  return(ret)
}

setMethod("[","MSMSpepXML",function(x,i,j="missing",drop="missing") "[.MSMSpepXML"(x,i))

"[.MSMSpepXML" <- function(x,i) {
  if (max(i)>nSpectrumQueries(x) | min(i)<1)
    stop("subscript out of bonds")
  if (length(i)==1)
    return(spectrumQueries(x)[[i]])
  return(spectrumQueries(x)[i])

}

setMethod("precNeutralMass","MSMSpepXML",
          function(object) return(sapply(spectrumQueries(object),precNeutralMass)))

setMethod("pepSequence","MSMSpepXML",
          function(object) lapply(spectrumQueries(object),pepSequence))

setMethod("scores","MSMSpepXML",
          function(object) lapply(spectrumQueries(object),scores))
