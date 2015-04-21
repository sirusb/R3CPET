## ----style, eval=TRUE, echo=FALSE, results='asis'---------------------------------------
BiocStyle::latex()

## ----options,echo=FALSE-----------------------------------
options(width=60)

## ----env, echo=FALSE, warning=FALSE-----------------------
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("igraph"))

## ----init-ChiapetExperimentData, results = 'hide', warning=FALSE, message=FALSE, prompt=TRUE, tidy=TRUE----
library(R3CPET)
petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  
x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 

## ----Read_chiapetData_chiapetTools, warning=FALSE, tidy=TRUE, prompt=TRUE----
petPath <- system.file("example","HepG2_interactions.txt",package="R3CPET")
petFile <- read.table(petPath, 
                      sep = "\t", header = TRUE)
head(petFile)

## ----Read_chiapetData_BED, warning=FALSE, tidy=TRUE, prompt=TRUE----
petPath <- system.file("example","HepG2_centered.bed",package="R3CPET")
petFile <- read.table(petPath, sep = "\t", header = FALSE, comment.char = '+')
                      
head(petFile)

## ----loadPETs_example, warning=FALSE, tidy=TRUE, prompt=TRUE----
## if it has 6 columns format IsBed = FALSE
petPath <- system.file("example","HepG2_interactions.txt",package="R3CPET")
x <- loadPETs(x,petFile=petPath, IsBed=FALSE)

## ----loadPETs_BED_example, warning=FALSE, tidy=TRUE, prompt=TRUE, eval=FALSE----
#  ## loading a 4 columns BED file
#  petPath <- system.file("example","HepG2_centered.bed",package="R3CPET")
#  x <- loadPETs(x,petFile=petPath, IsBed=TRUE, header = FALSE)

## ----TF_example, warning=FALSE, tidy=TRUE, prompt=TRUE----
## loading a 4 columns BED file
TFPath <- system.file("example","HepG2_TF.txt.gz",package="R3CPET")
TF <- read.table(TFPath, sep = "\t", header= FALSE)
head(TF)
x <- loadTFBS(x,tfbsFile=TFPath)

## ----HPRD_Biogrid_example, warning=FALSE, tidy=TRUE, prompt=TRUE----
data(HPRD)
data(Biogrid)

PPI.HPRD
PPI.Biogrid

## ----loadPPI_usage, warning=FALSE, tidy=TRUE, eval=FALSE, prompt=TRUE----
#  loadPPI(object,type=c("HPRD","Biogid"),customPPI = NULL,
#  	    filter = FALSE,term = "GO:0005634",
#  	     annot = NULL, RPKM = NULL, threshold = 1 )

## ----loadPPI_loading, warning=FALSE, tidy=TRUE, prompt=TRUE----
## loading the PPI with GO filtering
x <- loadPPI(x,type="HPRD", filter=TRUE)

## ----createIndexes, warning=FALSE, tidy=TRUE, prompt=TRUE----
x <- createIndexes(x)
x

## ----buildingNetworks, warning=FALSE, tidy=TRUE,cache=TRUE, prompt=TRUE, eval=FALSE----
#  nets<- buildNetworks(x, minFreq=0.1, maxFreq=0.9)
#  nets

## ----load_nets, warning=FALSE, echo=FALSE, message=TRUE,----
load(system.file("example","nets.RData",package="R3CPET"))
nets

## ----InferNetowrks_usage, warning=FALSE, tidy=TRUE, eval=FALSE----
#  InferNetworks(object,thr =0.5,max_iter = 500L, max_time = 3600L, ...)

## ----InferNetowrks, warning=FALSE, tidy=TRUE,cache=TRUE, prompt=TRUE, eval=FALSE----
#  hlda <- InferNetworks(nets)
#  hlda

## ----load_hlda, warning=FALSE, echo=FALSE, message=TRUE----
load(system.file("example","hlda.RData",package="R3CPET"))
hlda

## ----InferNetowrks_topElemets, warning=FALSE, tidy=TRUE, prompt=TRUE----
head(topEdges(hlda))
head(topNodes(hlda))

## ----InferNetowrks_Networks, warning=FALSE, tidy=TRUE, prompt=TRUE----
hlda <- GenerateNetworks(hlda)
head(networks(hlda))

## ----Annotate expression, warning=FALSE, tidy=TRUE, prompt=TRUE----
data(RPKMS)
hlda<- annotateExpression(hlda,RPKMS)
networks(hlda)[[1]]

## ----cluster_usage, tidy=TRUE,  eval=FALSE----------------
#  clusterInteractions(object, method=c("clues","sota"), nbClus=20 )

## ----clusterInteractions, tidy=TRUE,cache=TRUE, prompt=TRUE----
## clustering using the "clues" method
hlda <- clusterInteractions(hlda, method="clues")

## ----plot3CPETRes_usage, tidy=TRUE,  eval=FALSE-----------
#  
#  plot3CPETRes(object, path="", W=14, H=7 ,
#            type=c("heatmap","clusters","curve","avgCurve","netSim", "networks"),
#            byEdge=TRUE, layoutfct=layout.kamada.kawai, ...)

## ----getRegionsIncluster, tidy=TRUE, prompt=TRUE----------
getRegionsIncluster(hlda,x,cluster=3)

## ----heatmap, tidy=TRUE,cache=TRUE, prompt=TRUE,fig.align='center' , fig.height=3----
plot3CPETRes( hlda, type = 'heatmap')

## ----plot_curve, tidy=TRUE,cache=TRUE, prompt=TRUE,fig.align='center', fig.width=6, fig.height=3.5----
## plotting curves
plot3CPETRes(hlda, type = 'curve')

## ----plot_avgCurve, tidy=TRUE,cache=TRUE, prompt=TRUE,fig.align='center', fig.width=6, fig.height=3.5----
## plotting Average curves
plot3CPETRes(hlda, type = 'avgCurve')

## ----plot_clusters, tidy=TRUE,cache=TRUE, prompt=TRUE,fig.align='center',fig.width=5, fig.height=5----
## plotting pair-wise clusters scatter plots
plot3CPETRes(hlda, type = 'clusters')

## ----plot_networks, tidy=TRUE,cache=TRUE, prompt=TRUE,fig.align='center', warning=FALSE, message=FALSE, fig.width=3.5, fig.height=3.5----
nets_plot <- plot3CPETRes(hlda, type = 'networks')
plot(nets_plot[[4]])

## ----plot_netSim, tidy=TRUE,cache=FALSE, prompt=TRUE,fig.align='center', warning=FALSE, fig.width=5, fig.height=3.5----
plot3CPETRes(hlda,type = 'netSim')

## ----visualizeCircos_usage,tidy=TRUE,  eval=FALSE---------
#  
#  visualizeCircos(object, data, cluster = 1, chrLenghts = NULL)

## ----plot_circos, tidy=TRUE, prompt=TRUE,fig.align='center', fig.height=6, fig.width=6, message=FALSE, warning=FALSE----
visualizeCircos(hlda,x, cluster = 4)

## ----GOEnrich.networks_usage, warning=FALSE, tidy=TRUE, eval=FALSE, prompt=TRUE----
#  
#  GOEnrich.networks(object, fdr= 0.05, GOlimit = 5,path = "")

## ----GOEnrich.networks, warning=FALSE, tidy=TRUE, prompt=TRUE, fig.align='center', cache=TRUE, fig.show='hide', eval=FALSE----
#  
#  GOEnrich.networks(hlda, path= '.')

## ----outputGenesPerClusterToDir_usage,tidy=TRUE,  eval=FALSE----
#  outputGenesPerClusterToDir(hdaRes,data,path="ClustersGenes", ...)

## ----outputGenesPerClusterToDir, warning=FALSE, tidy=TRUE, prompt=TRUE, message=FALSE, eval=FALSE----
#  outputGenesPerClusterToDir(hlda, x)

## ----GOEnrich_folder, warning=FALSE, tidy=TRUE, prompt=TRUE, fig.align='center', fig.show='hide', message=FALSE, eval=FALSE----
#  GOEnrich.folder(folder = "ClustersGenes/")

## ----createShinyServer,tidy=TRUE,  eval=FALSE-------------
#  createServer(x,nets,hlda)

## ----sessioninfo------------------------------------------
sessionInfo()

