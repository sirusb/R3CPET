\name{GenerateNetworks-methods}
\docType{methods}
\alias{GenerateNetworks}
\alias{GenerateNetworks-methods}
\alias{GenerateNetworks,ChromMaintainers-method}
\title{ Generate a list of \code{igraph} networks}
\description{
 This methods converts the \code{networks} slot of a \code{\link{ChromMaintainers}} object, 
 it reads the \code{topEdge} slot and convert it into a list of \code{\link{igraph}} objectts.
}
\usage{
\S4method{GenerateNetworks}{ChromMaintainers}(object,\dots)
}
\arguments{
\item{object}{ a \code{\link{ChromMaintainers}} object}
\item{\dots}{future options not considered for the moment.}
}
\value{
 Returns \code{\link{ChromMaintainers}} object in which the \code{networks} slot is populated.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChromMaintainers}}, \code{\link{InferNetworks}}
}

\examples{

  ## get the different datasets path
  petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
  tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  
    
\dontrun{
  data(RPKMS)
  x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
  ## build the diffrent indexes
  x <- createIndexes(x)
  ## build networks connecting each interacting regions
  nets<- buildNetworks(x)

  ## infer the networks
  hlda<- InferNetworks(nets)

  hlda <- GenerateNetworks(hlda)
  networks(hlda)
}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
