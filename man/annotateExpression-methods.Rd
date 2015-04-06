\name{annotateExpression-methods}
\docType{methods}
\alias{annotateExpression}
\alias{annotateExpression-methods}
\alias{annotateExpression,ChromMaintainers,data.frame-method}

\title{Add the gene expression attribute to the graph nodes}
\description{
 This method is a kinda of helper method, it helps the user to add for each node in the inferred
 chromatin maintainer network the \code{RPKM} attributes. It is useful if the user want to save the 
 networks as \code{".gml"} files and visualize them using software such as Gephi or Cytoscape.
 Or maybe if he wants to know which networks are highly expressed then others.
}
\usage{
\S4method{annotateExpression}{ChromMaintainers,data.frame}(object, RPKMS)
}
\arguments{
\item{object}{a \code{ChromMaintainers} object in which the networks are already present}
\item{RPKMS}{ a two columns \code{data.frame}, where the first column contains the name of the gene 
and the second contains the expression values}
}
\value{
A \code{\link{ChromMaintainers}} object in which the networks are annotated.
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
  x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
  ## build the diffrent indexes
  x <- createIndexes(x)
  x
  
## build the different indexes
x <- createIndexes(x)

## build networks connecting each interacting regions
nets<- buildNetworks(x)

## infer the networks
hlda<- InferNetworks(nets)

networks(hlda)

## Annotate networks
hlda<- annotateExpression(hlda,as.data.frame(RPKMS))

## Notice the addition of the RPKM attribute to each network
networks(hlda)
}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
