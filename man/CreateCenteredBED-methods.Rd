\name{CreateCenteredBED-methods}
\docType{methods}
\alias{CreateCenteredBED-methods}
\alias{CreateCenteredBED,character-method}
\alias{CreateCenteredBED}
\title{ Create centered interactions}
\description{
 This helper method can be used to create a "bed" file in which the coordinates of the regions are the centre of the interactions
in the raw data. in  \pkg{R3CPET} we suppose that the centre of the interactions are the most enriched when doing read mapping
,thus, we consider just the region around the centre to detect the TFBS.
}
\usage{
\S4method{CreateCenteredBED}{character}(file, header=TRUE,dist=1000)
}
\arguments{
\item{file}{
a \code{character} indicating the location of the rawdata file. the file should be a six column "bed" file
in which the fist 3 columns indicate the left side interaction (chr, start, stop)
and the other 3 columns indicate the right side interaction.
}
\item{header}{
\code{logical}, indicates if the first line in the file is a header.
}
\item{dist}{
\code{numeric}, indicated the distance around the center of the region to take.
}
}
\value{
A 4 columns \code{data.frame} object, in which the first 3 columns indicate the location of the region and the 4th on indicate its name. 
The names are of the format \code{PET#\w+.1} for the left side regions and \code{PET#\w+.2} for the right side ones.
}
\references{
Mohamed Nadhir D, Yang C et al \emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, ....
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChiapetExperimentData}}, \code{\link{loadTFBS}} , \code{\link{loadPETs}}, \code{\link{loadPPI}}, \code{\link{createIndexes}}.
}
\examples{

## get interactions file location
petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")

res <- CreateCenteredBED(petFile, header=TRUE, dist=1000)
head(res)
}
\keyword{methods}
