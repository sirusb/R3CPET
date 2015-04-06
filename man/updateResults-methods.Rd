\name{updateResults-methods}
\docType{methods}
\alias{updateResults-methods}
\alias{updateResults,ChromMaintainers,NetworkCollection,numeric-method}
\alias{updateResults}
\title{update top maintain networks elements }
\description{
 This helper method can be used to update the list of interactions constituting the chromatin maintainer networks by changing the 
 threshold.
}
\usage{
\S4method{updateResults}{ChromMaintainers,NetworkCollection,numeric}(object,nets,thr=0.5)
}
\arguments{
\item{object}{a \code{\link{ChromMaintainers}} object in which the clusters are already calculated}
\item{nets}{a \code{\link{NetworkCollection}} used to get the vocabulary list.}
\item{thr}{The cut-off threshold. the top edges that have an enrichment value that sum up to a value >= \code{thr} are considered.}
}
\value{
 a \code{\link{ChromMaintainers}} object in which the \code{topNodes} and \code{topEdges} tables are updated.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{InferNetworks}}, \code{\link{NetworkCollection}}, \code{\link{ChromMaintainers}}
}
\examples{

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
    topNodes(hlda)

    hlda <- updateResults(hlda, nets, 0.4)
    topNodes(hlda)
    }
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
