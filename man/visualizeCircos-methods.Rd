\name{visualizeCircos-methods}
\docType{methods}
\alias{visualizeCircos}
\alias{visualizeCircos-methods}
\alias{visualizeCircos,ChromMaintainers,ChiapetExperimentData,numeric-method}
\title{ Generate circos plot per cluster}
\description{
 This method generates a basic circos plot of the chromatin interaction in a given cluster.
}
\usage{
\S4method{visualizeCircos}{ChromMaintainers,ChiapetExperimentData,numeric}(object, data, cluster = 1, chrLenghts = NULL)
}
\arguments{
\item{object}{ a \code{\link{ChromMaintainers}} object in which the clusters are already calculated. }
\item{data}{ a  \code{\link{ChiapetExperimentData}} containing the interaction data.}
\item{cluster}{the number of the cluster to display}
\item{chrLenghts}{the chromatin lengths. if not provided the package suppose it is a human chromatin and uses the 
corresponding lengths. Change it if you are using another species.}
}
\value{
\describe{
\item{\code{circos}}{ a \code{\link{GRanges}} object that contains the coordinate of the left side interactions. 
The right side interactions can be accessed by writing \code{circos\$to.gr}. }

\item{\code{plot}}{a \code{\link{ggplot}} object}
}
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{NetworkCollection}}, \code{\link{ChromMaintainers}} 
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
    hlda<- clusterInteractions(hlda)

    visualizeCircos(hlda,x, cluster=3)
}
}
\keyword{methods}
