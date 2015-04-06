\name{getRegionsInNetwork-methods}
\docType{methods}
\alias{getRegionsInNetwork-methods}
\alias{getRegionsInNetwork,ChromMaintainers,ChiapetExperimentData,numeric-method}
\alias{getRegionsInNetwork}
\title{ list of interactions per network}
\description{
 This method can be used to retrieve the genomic coordinated of the DNA-interactions
 enriched for each network given a certain threshold.
}
\usage{
\S4method{getRegionsInNetwork}{ChromMaintainers,ChiapetExperimentData,numeric}(hdaRes,data, net=1, thr=0.5, \dots)
}
\arguments{
\item{hdaRes}{a \code{\link{ChromMaintainers}} object which already contains the calculated results}
\item{data}{a \code{\link{ChiapetExperimentData}} object that contains the genomic location 
 of all the interactions.}
\item{net}{The ID of the network for which we want to get the list of the involved regions.}
\item{thr}{Specify the minimum enrichment value to consider an interaction to be controlled by the network. by default it is 0.5}
\item{\dots}{additional parameters not used for the moment.}
}
\value{
 a \code{\link{GRanges}} object is returned
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{InferNetworks}}, \code{\link{ChiapetExperimentData}}, \code{\link{ChromMaintainers}}
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

    ## infer the networks and do the clustering
    hlda<- InferNetworks(nets)

    ## return the DNA-interactions in cluster 3
    getRegionsIncluster(hlda,x,net=3)

}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }