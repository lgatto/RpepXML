## =================================================
## Class definition for a simplified element
## 'spectrum_hit' and children nodes of the pepXML standard.
## See pepXML specifications for underlying xsd
## -------------------------------------------------
setClass("SearchHit",
         representation(hitRank="integer",
                        proteins="character",
                        pepSequence="character",
                        modifications="data.frame",
                        scores="numeric"),
         contains = c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(searchHit="0.0.1")))
         )


## ===================
## Method definitions
## -------------------
setMethod("hitRank","SearchHit",function(object) object@hitRank)
setMethod("hitRank<-","SearchHit",
          function(object,value="integer") object@hitRank <- value)
setReplaceMethod("hitRank",
                 signature(object="SearchHit",
                           value="integer"),
                 function (object,value) {
                   object@hitRank <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("proteins","SearchHit",function(object) object@proteins)
setMethod("proteins<-","SearchHit",
          function(object,value="character") object@proteins <- value)
setReplaceMethod("proteins",
                 signature(object="SearchHit",
                           value="character"),
                 function (object,value) {
                   object@proteins <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("pepSequence","SearchHit",function(object) object@pepSequence)
setMethod("pepSequence<-","SearchHit",
          function(object,value="character") object@pepSequence <- value)
setReplaceMethod("pepSequence",
                 signature(object="SearchHit",
                           value="character"),
                 function (object,value) {
                   object@pepSequence <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("modifications","SearchHit",function(object) object@modifications)
setMethod("modifications<-","SearchHit",
          function(object,value="character") object@modifications <- value)
setReplaceMethod("modifications",
                 signature(object="SearchHit",
                           value="data.frame"),
                 function (object,value) {
                   object@modifications <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("scores","SearchHit",function(object) object@scores)
setMethod("scores<-","SearchHit",
          function(object,value="numeric") object@scores <- value)
setReplaceMethod("scores",
                 signature(object="SearchHit",
                           value="numeric"),
                 function (object,value) {
                   object@scores <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("ionScore","SearchHit",
          function(object) object@scores["ionscore"])

## =========================
## Other method definitions
## -------------------------
setMethod("show","SearchHit",
          function(object) {
            cat("Object of class",class(object),"\n")
            cat(" Rank:",hitRank(object),"\n")
            cat(" Sequence:",pepSequence(object),"\n")
            cat(" Scores:",scores(object),"\n")
            cat("",length(proteins(object)),"proteins:",
                paste(proteins(object),collapse=","),"\n")
            cat("",nrow(modifications(object)),"modifications\n")
          })



setMethod("filterHits","SearchHit",
          function(object,rank=1,ionscore=20) {
            if ((ionScore(object)>=ionscore) &&
                hitRank(object)<=rank)
              return(TRUE)
            return(FALSE)
          })

setMethod("nProteins","SearchHit",
          function(object) length(proteins(object))
          )
