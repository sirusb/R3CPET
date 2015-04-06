\name{outputGenesPerClusterToDir-methods}
\docType{methods}
\alias{outputGenesPerClusterToDir}
\alias{outputGenesPerClusterToDir-methods}
\alias{outputGenesPerClusterToDir,ChromMaintainers,ChiapetExperimentData-method}
\title{ List of genes in each cluster}
\description{
 This helper methods get the set of genes located in the DNA-regions in each cluster.
 A folder that contain a bunch of \code{.txt} files (one for each cluster) is generated.
 We consider a gene to be part of a cluster if the (-1500bp, +500bp) around its \code{TSS} 
 intersects with one of the DNA regions of the cluster.
 
}
\usage{
\S4method{outputGenesPerClusterToDir}{ChromMaintainers,ChiapetExperimentData}(hdaRes,data,path="ClustersGenes", \dots)
}
\arguments{
\item{hdaRes}{a \code{\link{ChromMaintainers}} object in which the DNA-interaction are already clustered.}
\item{data}{ a \code{\link{ChiapetExperimentData}} object containing information about the interactions.}
\item{path}{path of the folder to create. by default a folder named \code{ClustersGenes} is created
 in the current working directory.}
\item{\dots}{additional parameters, not used for the moment.}
}
\value{
The specified folder is created with a list \code{.txt} files that contain the list of genes.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChromMaintainers}}, \code{\link{InferNetworks}}, \code{\link{ChiapetExperimentData}}, 
\code{\link{clusterInteractions}}
}
\examples{

    petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
    tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  
    
\dontrun{
    x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
    ## build the different indexes
    x <- createIndexes(x)

    ## build networks connecting each interacting regions
    nets<- buildNetworks(x)

    ## infer the networks
    hlda<- InferNetworks(nets)
    hlda<- clusterInteractions(hlda)

    ## get the list of genes per cluster.
    outputGenesPerClusterToDir(hlda,x)

}
}
\keyword{methods}
