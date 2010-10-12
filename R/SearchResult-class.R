## =================================================
## Class definition for a simplified element
## 'search_result' of the pepXML standard.
## See pepXML specifications for underlying xsd
## -------------------------------------------------
setClass("SearchResult",
         representation(searchHits="list",
                        searchId="integer"),
         contains = c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(searchResults="0.0.1")))
         )

## ============================
## Accessor method definitions
## ----------------------------
setMethod("searchHits","SearchResult",function(object) object@searchHits)
setMethod("searchHits<-","SearchResult",
          function(object,value="list") object@searchHits <- value)
setReplaceMethod("searchHits",
                 signature(object="SearchResult",
                           value="list"),
                 function (object,value) {
                   object@searchHits <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("searchId","SearchResult",function(object) object@searchId)
setMethod("searchId<-","SearchResult",
          function(object,value="integer") object@searchId <- value)
setReplaceMethod("searchId",
                 signature(object="SearchResult",
                           value="integer"),
                 function (object,value) {
                   object@searchId <- value
                   if (validObject(object))
                     return(object)
                 })

## =========================
## Other method definitions
## -------------------------
setMethod("show","SearchResult",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat(" Search id:",searchId(object),"\n")
            cat("",nSearchHits(object),"search hit(s)\n")
          })

## =============================
## Accessors for SearchHit slots
## -----------------------------
setMethod("scores","SearchResult",
          function(object) {
            sapply(searchHits(object),
                   scores
                   )}
          )

setMethod("hitRank","SearchResult",
          function(object) {
            sapply(searchHits(object),
                   hitRank
                   )}
          )

setMethod("proteins","SearchResult",
          function(object) {
            unlist(sapply(searchHits(object),
                          proteins
                          ))}
          )

setMethod("pepSequence","SearchResult",
          function(object) {
            unlist(sapply(searchHits(object),
                          pepSequence
                          ))}
          )
## =========================
## Other method definitions
## -------------------------
setMethod("nSearchHits","SearchResult",
          function(object) length(searchHits(object))
          )

setMethod("filterHits","SearchResult",
          function(object,rank=1,ionscore=20) {
            sel <- sapply(searchHits(object),filterHits,rank,ionscore)
            if (sum(sel)>0) {
              searchHits(object) <- searchHits(object)[sel]
            } else {
              searchHits(object) <- list()
            }
            if (validObject(object))
              return(object)
          })

setMethod("nProteins","SearchResult",
          function(object) {
            if (nSearchHits(object)==0)
              return(0)
            sapply(searchHits(object),nProteins)
          })

setMethod("[","SearchResult",function(x,i,j="missing",drop="missing") "[.SearchResult"(x,i))

"[.SearchResult" <- function(x,i) {
  if (max(i)>nSearchHits(x) | min(i)<1)
    stop("subscript out of bonds")
  if (length(i)==1)
    return(searchHits(x)[[i]])
  return(searchHits(x)[i])
}
