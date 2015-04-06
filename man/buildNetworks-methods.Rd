\name{buildNetworks}
\docType{methods}
\alias{buildNetworks,ChiapetExperimentData-method}
\alias{buildNetworks}
\title{ Building interaction networks connecting interacting regions }
\description{
 This methods uses the background PPI to try to build an interaction network that connects 
 each interacting regions. If a \code{regionA} interacts with a \code{regionB} and if 
 \eqn{TF_A} is the list of TF in \code{regionA} and \eqn{TF_B} is the list of TF in 
 \code{regionB}, than we use the loaded PPI as a background network to connect each TF from 
 \eqn{TF_A} to each TF in \eqn{TF_B}. 
 
 We suppose that a minimum number of physical interactions (minimum energy) are needed to connection
 each TF to the other. Thus, we take the shortest path in the PPI. at this stage, each network is a collection
 of edges. 
}
\usage{
\S4method{buildNetworks}{ChiapetExperimentData}(object, minFreq = 0.25, maxFreq = 0.75)
}
\arguments{
\item{object}{
 a \code{\link{ChiapetExperimentData}} object in which the interactions and TFBS and PPI are already loaded. 
 Check \code{\link{loadPETs}}, \code{\link{loadTFBS}}, \code{\link{loadPPI}} for more info.
}
\item{minFreq}{
 After constructing the networks for all the interacting regions all edges that appear in less than \code{minFreq} of the networks are considered to be outliers.
}
\item{maxFreq}{
 After constructing the networks for all the interacting regions all edges that appear in more than \code{maxFreq} of the networks are considered to be interactions involving general TF and are removed.
}
}
\value{
A \code{\link{NetworkCollection}} object that contain the list of all the constructed networks and their sizes.

NOTE: interactions for which no TF was bound or no networks could be constructed or which was empty 
after filtering will not be considered.
}
\references{
Mohamed Nadhir D, Yang C et al \emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, ....
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChiapetExperimentData}}, \code{\link{loadTFBS}} , \code{\link{loadPETs}}, \code{\link{loadPPI}}, \code{\link{createIndexes}}
}
\examples{
    
  ## get the different datasets path
  petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
  tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  
    
  
   \dontrun{  
  x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
  ## build the diffrent indexes
  x <- createIndexes(x)  

  ## build networks connecting each interacting regions
  nets<- buildNetworks(x)
  nets
}
}
\keyword{methods}
