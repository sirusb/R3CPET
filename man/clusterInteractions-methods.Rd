\name{clusterInteractions-methods}
\docType{methods}
\alias{clusterInteractions}
\alias{clusterInteractions,ChromMaintainers-method}
\title{ Grouping DNA interactions by enrichment profile }
\description{
 This method aims at clustering the DNA interactions according to their partnership probability to 
 the inferred chromatin maintainer networks.
 
 Two kinds of clustering are supported supervised and un-supervised. In the first one the function\code{\link{sota}}
 from the \code{\link{clValid}} package. In the second case the \code{\link{clues}} method from the \code{\link{clues}}
 package is used.
}
\usage{
\S4method{clusterInteractions}{ChromMaintainers}(object, method=c("clues","sota"), nbClus=20 )
}
\arguments{
 \item{object}{ (\emph{\strong{Required})} a non-empty \code{\link{ChromMaintainers}} object}
 \item{method}{
 (\emph{optional})used to specify the method to use. if \code{clues} is specified, the number of clusters will be 
 determined automatically. if \code{method = "sota"} , the user can specify the number of clusters
 by setting the parameter \code{nbClus}, by default it is set to 20. by default \code{method = "clues"}.
 
 Another option is that user can first do an automatic clustering using \code{"clues"} and then if judges
 that more or less clusters are required he can recall the \code{clusterInteractions} method and specify the number the 
 number of clusters.
 }
 \item{nbClus}{
 (\emph{optional}) The user-specified number of clusters. It is taken into consideration only if \code{method = sota}.
 }
}
\value{
A \code{\link{ChromMaintainers}} object in which the \code{clusRes} is populated as a \code{\link{sota}} 
or \code{\link{clues}} object.
}
\references{
Herrero, J., Valencia, A, and Dopazo, J. (2005).  \emph{A hierarchical unsupervised growing neural network for clustering gene expression patterns}. Bioinformatics, 17, 126-136.

Wang, S., Qiu, W., and Zamar, R. H. (2007). \emph{ CLUES: A non-parametric clustering method based on local shrinking}. Computational Statistics & Data Analysis, Vol. 52, issue 1, pages 286-298
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChromMaintainers}}, \code{\link{clues}} , \code{\link{sota}}, \code{\link{InferNetworks}}
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
