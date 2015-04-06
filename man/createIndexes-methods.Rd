\name{createIndexes-methods}
\docType{methods}
\alias{createIndexes-methods}
\alias{createIndexes,ChiapetExperimentData-method}
\alias{createIndexes}
\title{ Preparing TF indexes per region}
\description{
 After loading the interactions and the TFBS, the \code{createIndexes} method can be used to build indexes for fast 
 look-up for which which TF are located in which region. This method is an intermediate step needed for further analysis.
}
\usage{
\S4method{createIndexes}{ChiapetExperimentData}(object, minOverlap = 50)
}
\arguments{
 \item{object}{
 a \code{\link{ChiapetExperimentData}} object in which the interactions and TFBS are already loaded. 
 Check \code{\link{loadPETs}} and \code{\link{loadTFBS}} for more info.
 }
  \item{minOverlap}{
 The minimum overlap between a TF binding site and a region, to consider a TF as binding to that region. 
 The default value is \code{50}.
 }
}
\value{
A \code{\link{ChiapetExperimentData}} object in which the \code{.dt} slot is populated as a \code{\linkS4class{data.table}} object.
}
\references{
Mohamed Nadhir D, Yang C et al \emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, ....
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChiapetExperimentData}}, \code{\link{loadTFBS}} , \code{\link{loadPETs}}, \code{\link{loadPPI}}
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
}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
