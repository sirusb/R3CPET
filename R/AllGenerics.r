setGeneric("loadPETs",signature=c("object","petFile"),
           function(object, petFile,
                    IsBed=TRUE, header=TRUE,
                    dist=1000) standardGeneric("loadPETs"))

setGeneric("pet", signature = c( "object"),
           function(object) standardGeneric("pet"))

setGeneric("pet<-", function(object, value) standardGeneric("pet<-"))
               
setGeneric("tfbs", signature = c("object"),
           function(object) standardGeneric("tfbs"))

setGeneric("tfbs<-", function(object, value) standardGeneric("tfbs<-"))
    
setGeneric("loadTFBS", function(object,tfbsFile, ...)
  standardGeneric("loadTFBS"))

setGeneric("visualizeInteractions", signature=c("object","range"),
           function(object, range) standardGeneric("visualizeInteractions"))

setGeneric("plotTrack", signature=c("object","range"),
           function(object, range) standardGeneric("plotTrack"))

setGeneric("loadPPI",signature=c("object"),
           function(object,type=c("HPRD","Biogid"),
                    customPPI= NULL, 
                    filter = FALSE,
                    term ="GO:0005634",
                    annot=NULL,
                    RPKM= NULL, threshold=1 ) standardGeneric("loadPPI"))

setGeneric("ppi", signature = c("object"),
           function(object) standardGeneric("ppi"))

setGeneric("ppi<-", function(object, value) standardGeneric("ppi<-"))
               
setGeneric("createIndexes", signature=c("object"),
           function(object, ...) standardGeneric("createIndexes"))
          
          
setGeneric("PrepareData", function(petFile,tfbsFile, petIsBed)
           standardGeneric("PrepareData"))


setGeneric("CreateCenteredBED",signature =c("file"),
           function(file, ...)  standardGeneric("CreateCenteredBED"))


setGeneric("buildNetworks", signature= c("object"),
          function(object,...)
                   #minFreq,
                   #maxFreq)
                   standardGeneric("buildNetworks") )

setGeneric("networks", signature=  c("object"),
           function(object) standardGeneric("networks"))

setGeneric("sizes", signature = c("object"),
           function(object) standardGeneric("sizes"))

setGeneric("TF", signature = c("object"),
           function(object) standardGeneric("TF"))

setGeneric("InferNetworks", signature= c("object"),
           function(object,...) standardGeneric("InferNetworks"))

setGeneric("clusterInteractions", signature= c("object"),
          function(object, ...) standardGeneric("clusterInteractions"))

setGeneric("getClusters", signature= c("object"),
           function(object) standardGeneric("getClusters"))

setGeneric("getRegionsIncluster", signature=c("hdaRes", "data","cluster"),
          function(hdaRes, data,cluster, ...) standardGeneric("getRegionsIncluster"))

setGeneric("getRegionsInNetwork", signature=c("hdaRes", "data","net"),
           function(hdaRes, data,net, ...) standardGeneric("getRegionsInNetwork"))

setGeneric("plot3CPETRes", signature= c("object"),
           function(object, ...) standardGeneric("plot3CPETRes"))

setGeneric("topEdges", signature= c("object"),
           function(object)  standardGeneric("topEdges"))

setGeneric("topNodes", signature = c("object"),
           function(object) standardGenetic("topNodes"))

setGeneric("wordsPerTopic", signature = c("object"),
           function(object) standardGneric("wordsPerTopic"))

setGeneric("docPerTopic", signature = c("object"),
           function(object) standrdGeneric("docPerTopic"))

setGeneric("betas", signature = c("object"),
           function(object) standardGeneric("betas"))

setGeneric("updateResults", signature=c("object","nets", "thr"),
           function(object, nets, thr) standardGeneric("updateResults"))

setGeneric("outputGenesPerClusterToDir", signature= c("hdaRes", "data"),
           function(hdaRes, data, ...) standardGeneric("outputGenesPerClusterToDir"))

setGeneric("outputGenesPerNetworkToDir", signature= c("hdaRes", "data"),
           function(hdaRes, data, ...) standardGeneric("outputGenesPerNetworkToDir"))

setGeneric("GOEnrich.folder",signature = c("folder"),
           function(folder, ...) standardGeneric("GOEnrich.folder"))

setGeneric("GOEnrich.networks",signature= c("object"),
           function(object, ...) standardGeneric("GOEnrich.networks"))

setGeneric("GenerateNetworks", signature = c("object"),
           function(object, ...) standardGeneric("GenerateNetworks") )

setGeneric("annotateExpression", signature= c("object", "RPKMS"),
           function(object, RPKMS) standardGeneric("annotateExpression"))


setGeneric("visualizeCircos", signature = c("object", "data",  "cluster"),
           function(object, data, cluster, ...) standardGeneric("visualizeCircos"))


setGeneric("createServer", signature = c("x","nets","hlda"),
          function(x, nets, hlda) standardGeneric("createServer"))
