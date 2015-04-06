\name{PrepareData-methods}
\docType{methods}
\alias{PrepareData-methods}
\alias{PrepareData,character,character,logical-method}
\alias{PrepareData}
\title{Loading the raw data all at once}
\description{
 Instead of loading the data one at a time and then creating the index using the methods \code{\link{loadPETs}},\code{\link{loadTFBS}} and \code{\link{createIndexes}}. The user can directly use the method \code{PrepareData} to do that.
}
\usage{
\S4method{PrepareData}{character,character,logical}(petFile,tfbsFile, petIsBed=TRUE)
}
\arguments{
  \item{petFile}{
  a \code{character} specifying the path to the interaction file. it the file is in a "bed" format \code{petIsBed} should be \code{TRUE}. 
  The data should be formated as described in \code{\link{loadPETs}}.
  }
  \item{tfbsFile}{
  a \code{character} specifying the path to the transcription factors binding site file. The data should be formated as described in \code{\link{loadTFBS}}.
  }
  \item{petIsBed}{
  a logical value specifying if the interaction file is in a "bed" format or not.
  }
}
\value{
A \code{\link{ChiapetExperimentData}} object in which the \code{pet},\code{tfbs} and \code{.dt} slots populated .
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

    ## get the TFBS file location
    tfFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")

    \dontrun{
    ## load the data
    x<- PrepareData(petFile, tfFile, FALSE)
    x
    }
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
