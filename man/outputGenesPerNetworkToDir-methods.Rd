\name{outputGenesPerNetworkToDir}
\docType{methods}
\alias{outputGenesPerNetworkToDir}
\alias{outputGenesPerNetworkToDir-methods}
\alias{outputGenesPerNetworkToDir,ChromMaintainers,ChiapetExperimentData-method}
\title{ List of genes controlled by each network}
\description{
 This helper methods get the set of genes located in the DNA-regions controlled by each network.
 A folder that contains a bunch of \code{.txt} files (one for each network) is generated.
 We consider (-2500bp, +2500bp) around the \code{TSS} of gene located in a region showing 0.5 or more enrichment 
 for the network.
}
\usage{
\S4method{outputGenesPerNetworkToDir}{ChromMaintainers,ChiapetExperimentData}(hdaRes,data,path="NetworksGenes", \dots)
}
\arguments{
\item{hdaRes}{a \code{\link{ChromMaintainers}} object containing the enrichment values.}
\item{data}{ a \code{\link{ChiapetExperimentData}} object containing information about the interactions.}
\item{path}{path of the folder to create. by default a folder named \code{NetworksGenes} is created
 in the current working directory.}
\item{\dots}{additional parameters, not used for the moment.}
}
\value{
The specified folder is created with a list \code{.txt} files each for each network that contain the list of genes.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChromMaintainers}}, \code{\link{InferNetworks}}, \code{\link{ChiapetExperimentData}}
}
\examples{
  ## get the different datasets path
  petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
  tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")    
\dontrun{
    data(RPKMS)
    x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
    ## build the different indexes
    x <- createIndexes(x)

    ## build networks connecting each interacting regions
    nets<- buildNetworks(x)

    ## infer the networks
    hlda<- InferNetworks(nets)

    ## get the list of genes per network.
    outputGenesPerNetworkToDir(hlda,x)

}
}
\keyword{methods}
