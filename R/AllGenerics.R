## ===============================
## Accessors for MSMSpepXML slots
## -------------------------------
setGeneric("pepFile",function(object) standardGeneric("pepFile"))
setGeneric("pepFile<-",function(object,value) standardGeneric("pepFile<-"))
setGeneric("searchEngine",function(object) standardGeneric("searchEngine"))
setGeneric("searchEngine<-",function(object,value) standardGeneric("searchEngine<-"))
setGeneric("sampleEnzyme",function(object) standardGeneric("sampleEnzyme"))
setGeneric("sampleEnzyme<-",function(object,value) standardGeneric("sampleEnzyme<-"))
setGeneric("searchDatabase",function(object) standardGeneric("searchDatabase"))
setGeneric("searchDatabase<-",function(object,value) standardGeneric("searchDatabase<-"))
setGeneric("spectrumQueries",function(object) standardGeneric("spectrumQueries"))
setGeneric("spectrumQueries<-",function(object,value) standardGeneric("spectrumQueries<-"))

## ==============================
## Accessors for SearchHit slots 
## ------------------------------
setGeneric("hitRank",function(object) standardGeneric("hitRank"))
setGeneric("hitRank<-",function(object,value) standardGeneric("hitRank<-"))
setGeneric("proteins",function(object) standardGeneric("proteins"))
setGeneric("proteins<-",function(object,value) standardGeneric("proteins<-"))
setGeneric("pepSequence",function(object) standardGeneric("pepSequence"))
setGeneric("pepSequence<-",function(object,value) standardGeneric("pepSequence<-"))
setGeneric("modifications",function(object) standardGeneric("modifications"))
setGeneric("modifications<-",function(object,value) standardGeneric("modifications<-"))
setGeneric("scores",function(object) standardGeneric("scores"))
setGeneric("scores<-",function(object,value) standardGeneric("scores<-"))
setGeneric("ionScore",function(object) standardGeneric("ionScore"))

## =================================
## Accessors for SearchResult slots
## ---------------------------------
setGeneric("searchHits",function(object) standardGeneric("searchHits"))
setGeneric("searchHits<-",function(object,value) standardGeneric("searchHits<-"))
setGeneric("searchId",function(object) standardGeneric("searchId"))
setGeneric("searchId<-",function(object,value) standardGeneric("searchId<-"))

## ==================================
## Accessors for SpectrumQuery slots
## ----------------------------------
setGeneric("searchResults",function(object) standardGeneric("searchResults"))
setGeneric("searchResults<-",function(object,value) standardGeneric("searchResults<-"))
setGeneric("spectrumId",function(object) standardGeneric("spectrumId"))
setGeneric("spectrumId<-",function(object,value) standardGeneric("spectrumId<-"))
setGeneric("startScan",function(object) standardGeneric("startScan"))
setGeneric("startScan<-",function(object,value) standardGeneric("startScan<-"))
setGeneric("endScan",function(object) standardGeneric("endScan"))
setGeneric("endScan<-",function(object,value) standardGeneric("endScan<-"))
setGeneric("precNeutralMass",function(object) standardGeneric("precNeutralMass"))
setGeneric("precNeutralMass<-",function(object,value) standardGeneric("precNeutralMass<-"))
setGeneric("assumedCharge",function(object) standardGeneric("assumedCharge"))
setGeneric("assumedCharge<-",function(object,value) standardGeneric("assumedCharge<-"))
setGeneric("queryIndex",function(object) standardGeneric("queryIndex"))
setGeneric("queryIndex<-",function(object,value) standardGeneric("queryIndex<-"))

## ==========================
## Other generic definitions
## --------------------------
setGeneric("nProteins",function(object) standardGeneric("nProteins"))
setGeneric("nSearchResults",function(object) standardGeneric("nSearchResults"))
setGeneric("nSearchHits",function(object) standardGeneric("nSearchHits"))
setGeneric("nSpectrumQueries",function(object) standardGeneric("nSpectrumQueries"))
setGeneric("filterHits",function(object,rank=NULL,ionscore=NULL) standardGeneric("filterHits"))
