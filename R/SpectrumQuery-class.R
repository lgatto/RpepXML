## =================================================
## Class definition for a simplified element
## 'spectrum_query' of the pepXML standard.
## See pepXML specifications for underlying xsd
## -------------------------------------------------
setClass("SpectrumQuery",
         representation(searchResults="list",
                        spectrumId="character",
                        startScan="integer",
                        endScan="integer",
                        precNeutralMass="numeric",
                        assumedCharge="numeric",
                        queryIndex="integer"
                        ),
         contains = c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(spectrumQuery="0.0.1")))
         )


## ============================
## Accessor method definitions
## ---------------------------
setMethod("searchResults","SpectrumQuery",function(object) object@searchResults)
setMethod("searchResults<-","SpectrumQuery",
          function(object,value="list") object@searchResults <- value)
setReplaceMethod("searchResults",
                 signature(object="SpectrumQuery",
                           value="list"),
                 function(object,value) {
                   object@searchResults <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("spectrumId","SpectrumQuery",function(object) object@spectrumId)
setMethod("spectrumId<-","SpectrumQuery",
          function(object,value="character") object@searchResults <- value)
setReplaceMethod("spectrumId",
                 signature(object="SpectrumQuery",
                           value="character"),
                 function(object,value) {
                   object@spectrumId <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("startScan","SpectrumQuery",function(object) object@startScan)
setMethod("startScan<-","SpectrumQuery",
          function(object,value="integer") object@searchResults <- value)
setReplaceMethod("startScan",
                 signature(object="SpectrumQuery",
                           value="integer"),
                 function(object,value) {
                   object@startScan <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("endScan","SpectrumQuery",function(object) object@endScan)
setMethod("endScan<-","SpectrumQuery",
          function(object,value="integer") object@searchResults <- value)
setReplaceMethod("endScan",
                 signature(object="SpectrumQuery",
                           value="integer"),
                 function(object,value) {
                   object@endScan <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("precNeutralMass","SpectrumQuery",function(object) object@precNeutralMass)
setMethod("precNeutralMass<-","SpectrumQuery",
          function(object,value="numeric") object@searchResults <- value)
setReplaceMethod("precNeutralMass",
                 signature(object="SpectrumQuery",
                           value="numeric"),
                 function(object,value) {
                   object@precNeutralMass <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("assumedCharge","SpectrumQuery",function(object) object@assumedCharge)
setMethod("assumedCharge<-","SpectrumQuery",
          function(object,value="numeric") object@searchResults <- value)
setReplaceMethod("assumedCharge",
                 signature(object="SpectrumQuery",
                           value="numeric"),
                 function(object,value) {
                   object@assumedCharge <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("queryIndex","SpectrumQuery",function(object) object@queryIndex)
setMethod("queryIndex<-","SpectrumQuery",
          function(object,value="numeric") object@queryIndex <- value)
setReplaceMethod("queryIndex",
                 signature(object="SpectrumQuery",
                           value="integer"),
                 function(object,value) {
                   object@queryIndex <- value
                   if (validObject(object))
                     return(object)
                 })

## =========================
## Other method definitions
## -------------------------
setMethod("show","SpectrumQuery",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat(" Spectrum id:",spectrumId(object),"\n")
            cat(" Start - end scans:",startScan(object),"-",endScan(object),"\n")
            cat(" Precursor neutral mass:",precNeutralMass(object),"\n")
            cat(" Assumed charge:",assumedCharge(object),"\n")
            cat(" Query index:",queryIndex(object),"\n")
            cat("",nSearchResults(object),"search result(s)\n")
            if (nSearchResults(object)>0) {
              cat("",sum(nSearchHits(object)),"search hits(s)\n")
              cat("",sum(unlist(nProteins(object))),"proteins(s)\n")
            }
          })

setMethod("nSearchResults","SpectrumQuery",
          function(object) length(searchResults(object))
          )

setMethod("nSearchHits","SpectrumQuery",
          function(object) {
            if (nSearchResults(object)>0) {
              return(sapply(searchResults(object),nSearchHits))
            } 
            return(0)
          })

setMethod("proteins","SpectrumQuery",
          function(object) {
            as.character(unlist(sapply(searchResults(object),proteins)))
          })

setMethod("filterHits","SpectrumQuery",
          function(object,rank=1,ionscore=20) {
            ## If there are no search results
            ## return object as is
            if (nSearchResults(object)==0)
              return(object)
            ## If there are search results,
            ## each one is filtered
            fsrchreslst <- lapply(searchResults(object),filterHits,rank,ionscore)
            sel <- sapply(fsrchreslst,function(x) {
              nSearchHits(x)>0
            })
            ## If some search results still have search hits,
            ## keep only those, else update search results
            ## as an empty list
            if (sum(sel)>0) {
              searchResults(object) <- fsrchreslst[sel]
            } else {
              searchResults(object) <- list()
            }
            if (validObject(object))
              return(object)
          })

setMethod("nProteins","SpectrumQuery",
          function(object) {
            if (length(searchResults(object))==0)
              return(0)
            return(sapply(searchResults(object),nProteins))
          })

setMethod("[","SpectrumQuery",function(x,i,j="missing",drop="missing") "[.SpectrumQuery"(x,i))

"[.SpectrumQuery" <- function(x,i) {
  if (max(i)>nSearchResults(x) | min(i)<1)
    stop("subscript out of bonds")
  if (length(i)==1)
    return(searchResults(x)[[i]])
  return(searchResults(x)[i])
}


setMethod("pepSequence","SpectrumQuery",
          function(object) lapply(searchResults(object),pepSequence))

setMethod("scores","SpectrumQuery",
          function(object) lapply(searchResults(object),scores))

