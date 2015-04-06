\name{createServer-methods}
\docType{methods}
\alias{createServer-methods}
\alias{createServer,ChiapetExperimentData,NetworkCollection,ChromMaintainers-method}
\alias{createServer}
\title{ Explore results in a web browser}
\description{
 To facilitate the interaction of the user with the package, we added an additional web interface using the \pkg{shiny} package.
 The user can check some basic statistics about the row data and he can explore the results and generate some graphs.
}
\usage{
\S4method{createServer}{ChiapetExperimentData,NetworkCollection,ChromMaintainers}(x,nets,hlda)
}
\arguments{
\item{x}{a \code{\link{ChiapetExperimentData}} object in which the interactions, TFBS and the index tables are already created.}
\item{nets}{a \code{\link{NetworkCollection}} object containing the list of the used TF and the initial interaction networks.}
\item{hlda}{a \code{\link{ChromMaintainers}} object in which the results are already calculated.}
}
\value{
 A webpage is opened.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{NetworkCollection}}, \code{\link{ChiapetExperimentData}}, \code{\link{ChromMaintainers}}
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

  ## infer the networks and do the clustering 
  hlda<- InferNetworks(nets)
  hlda<- clusterInteractions(hlda)

  ## Run the server
  createServer(x, nets, hlda)

}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
