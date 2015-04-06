\name{getRegionsIncluster-methods}
\docType{methods}
\alias{getRegionsIncluster}
\alias{getRegionsIncluster-methods}
\alias{getRegionsIncluster,ChromMaintainers,ChiapetExperimentData,numeric-method}
\title{ list of interactions per cluster }
\description{
 This method can be used to retrieve the genomic coordinated of the DNA-interactions
 in each cluster.
}
\usage{
\S4method{getRegionsIncluster}{ChromMaintainers,ChiapetExperimentData,numeric}(hdaRes,data, cluster=1, \dots)
}
\arguments{
\item{hdaRes}{a \code{\link{ChromMaintainers}} object in which the clusters are already calculated}
\item{data}{a \code{\link{ChiapetExperimentData}} object that contains the genomic location 
 of all the interactions.}
\item{cluster}{The ID of the cluster for which we want to get the list of the involved regions.}
\item{\dots}{additional parameters not used for the moment.}
}
\value{
 a \code{\link{GRanges}} object is returned
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{clusterInteractions}}, \code{\link{InferNetworks}},
 \code{\link{ChiapetExperimentData}}, \code{\link{ChromMaintainers}}
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
    hlda<- clusterInteractions(hlda)

    ## return the DNA-interactions in cluster 3
    getRegionsIncluster(hlda,x,cluster=3)

}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
