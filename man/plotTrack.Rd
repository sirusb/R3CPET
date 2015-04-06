\name{plotTrack-methods}
\docType{methods}
\alias{plotTrack}
\alias{plotTrack-methods}
\alias{plotTrack,ChiapetExperimentData,GRanges-method}
\title{ Plot interaction on a genomic track }
\description{
 This helper method can be used to display a genomic track for a certain location that contains the chromosome and the 
 related interactions if any.
}
\usage{
\S4method{plotTrack}{ChiapetExperimentData,GRanges}(object, range)
}
\arguments{
\item{object}{a \code{\link{ChiapetExperimentData}} object which contains the raw data.}
\item{range}{ The genomic coordinates of the location to display}
}
\value{
 a \code{ggbio::track} object
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChiapetExperimentData}}
}
\examples{

    petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
    tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  

    \dontrun{
    x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
    gr <- GRanges("chr1",IRanges(start=100000,end=300000))
    plotTrack(x,gr)
}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
