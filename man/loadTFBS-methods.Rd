\name{loadTFBS-methods}
\docType{methods}
\alias{loadTFBS-methods}
\alias{loadTFBS,ChiapetExperimentData}
\alias{loadTFBS,ChiapetExperimentData,character-method}
\alias{loadTFBS}

\title{ Loading TF binding sites}
\description{
This methods reads a BED file that contains the peak positions of different TF. 
All the TF peaks should be merged into one BED file that contains 4 columns that respectively contain
the chromosome name, peak start, peak end, TF name.
}
\usage{
\S4method{loadTFBS}{ChiapetExperimentData,character}(object, tfbsFile,header=FALSE, \dots)
}
\arguments{
\item{object}{
 (\strong{Required}) a \code{\link{ChiapetExperimentData}} object

 }
 \item{tfbsFile}{
 (\strong{Required}): path the BED file containing the position of the different TF binding site
 All the TF binding sites should be merged in this file as showed in this example:
 \preformatted{
 chr1  569820  569998 BHLHE40
 chr1  936071  936346 BHLHE40
 chr1 1014795 1015082 BHLHE40
 ...........................
 ...........................
chrY 13485240 13485769 ZBTB33
chrY 13488718 13489030 ZBTB33
chrY 15016340 15016848 ZBTB33
chrY 58843918 58844104 ZBTB33
 }
 }
 \item{header}{
 \emph{(optional)} indicates if the provided BED file has a header or node. by default \code{header=FALSE}
 },
 \item{\dots}{ reserved for later use.}
}
\value{
A \code{\link{ChiapetExperimentData}} object in which the \code{tfbs} slot is populated as a \code{\linkS4class{GRanges}} object.
}
\references{
Mohamed Nadhir D, Yang C et al \emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, ....
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChiapetExperimentData}}, \code{\link{loadTFBS}} , \code{\link{loadPPI}}, \code{\link{createIndexes}}
}
\examples{
## Create a ChiapetExperimentData object
x <- ChiapetExperimentData(ppiType= "HPRD")

## load TFBS
tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")
x <- loadTFBS(x,tfbsFile=tfbsFile)
tfbs(x)
}
\keyword{methods}
