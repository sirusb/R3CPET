\name{clusterInteractions-methods}
\docType{methods}
\alias{clusterInteractions}
\alias{clusterInteractions,ChromMaintainers-method}
\title{ Grouping DNA interactions by enrichment profile }
\description{
 This method aims at clustering the DNA interactions according to their partnership probability to 
 the inferred chromatin maintainer networks. 
}
\usage{
\S4method{clusterInteractions}{ChromMaintainers}(object, method="sota", nbClus=20 )
}
\arguments{
 \item{object}{ (\emph{\strong{Required})} a non-empty \code{\link{ChromMaintainers}} object}
 \item{method}{
 (\emph{optional})used to specify the method to use. Only the  \code{method = "sota"} is supported for the moment,
the \code{method='clues'} is deprecated. The user needs to specify the number of clusters
 by setting the parameter \code{nbClus}, by default it is set to 20.
 
 }
 \item{nbClus}{
 (\emph{optional}) The user-specified number of clusters. It is taken into consideration only if \code{method = sota}.
 }
}
\value{
A \code{\link{ChromMaintainers}} object in which the \code{clusRes} is populated as a \code{\link{sota}}.
}
\references{
Herrero, J., Valencia, A, and Dopazo, J. (2005).  \emph{A hierarchical unsupervised growing neural network for clustering gene expression patterns}. Bioinformatics, 17, 126-136.

}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChromMaintainers}},  \code{\link{sota}}, \code{\link{InferNetworks}}
}
\examples{

 data(RPKMS)

  ## get the different datasets path
  petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
  tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  
    
\dontrun{
  x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
  ## build the diffrent indexes
  x <- createIndexes(x)
  x

  ## build networks connecting each interacting regions
  nets<- buildNetworks(x)

  ## infer the networks
  hlda<- InferNetworks(nets)

  #cluster
  hlda<- clusterInteractions(hlda)

  #Display heatmap
  plot3CPETRes(hlda,type="heatmap")
  hlda
}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
