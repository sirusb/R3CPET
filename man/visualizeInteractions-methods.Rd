\name{visualizeInteractions}
\docType{methods}
\alias{visualizeInteractions-methods}
\alias{visualizeInteractions,ChiapetExperimentData,GRanges-method}
\alias{visualizeInteractions}
\title{ Display a Circos plot of ChIA-pet interactions  }

\description{
This method can be used to draw a circos plot of the chromatin interactions located 
in the given genomic range.
}
\usage{
\S4method{visualizeInteractions}{ChiapetExperimentData,GRanges}(object, range)
}
\arguments{
  \item{object}{ a \code{\link{ChiapetExperimentData}} in which the interactions are already loaded.
  Check \code{\link{loadPETs}} for more info.
}
  \item{range}{ a \code{\link{GRanges}} object containing the genomic region of interest.}
}
\value{
A ciros plot of the selected region is displayed and a list containing the following objects is returned.
\describe{
    \item{\code{circos} :}{a \code{\link{GRanges}} object that contains the involved chromatin interactions.}
    \item{\code{plot} :}{a \code{\link{ggplot}} object containing plot.}    
   }
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChiapetExperimentData}}, \code{\link{loadPETs}}, \code{\link{ggbio}}, \code{\link{GRanges}}
}
\examples{
  
  petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
  tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  
  
  
  \dontrun{
  x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
  ## plot intractions in the region of interest
  gr <- GRanges("chr1", IRanges(1240000,10300000))
  p <- visualizeInteractions(x, gr)
  p
  }
}
\keyword{methods}
