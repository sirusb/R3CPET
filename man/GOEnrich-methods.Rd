\name{GOEnrich.folder-methods}
\docType{methods}
\alias{GOEnrich.folder}
\alias{GOEnrich.folder-methods}
\alias{GOEnrich.folder,character-method}
\alias{GOEnrich.networks}
\alias{GOEnrich.networks-methods}
\alias{GOEnrich.networks,ChromMaintainers-method}

\title{GO enrichment methods}
\description{
 This helper methods can be called to do GO enrichment by using the \code{DAVID} web service.
 
 \code{GOEnrich.networks} can be used to do a GO enrichment of the chromatin maintainer networks.
 
 \code{GOEnrich.folder} can be called to do a GO enrichment on the gene-list files generated
 by the method \code{outputGenesPerClusterToDir}.  
 
 There is a \code{5 secs} delay between each request to not avoid being rejected 
 by the server.
}
\usage{
\S4method{GOEnrich.folder}{character}(folder, fdr=0.05,GOlimit=20)

\S4method{GOEnrich.networks}{ChromMaintainers}(object, fdr=0.05, GOlimit= 5,path="")
}
\arguments{
\item{folder}{ name of the folder that contains the gene-list files. The files are supposed to 
 have a \code{.txt} extension. The first column of each file is supposed to contain the genes EntrezID. }
 \item{object}{ a \code{"ChromMaintainers"} objects with the \code{topNodes} already calculated.}
\item{fdr}{ cut-off value GO terms with fdr value <= fdr will be considered. 
     Benjamini-Hochberg FDR is used.}
\item{GOlimit}{ the number of top GO terms to return.}
\item{path}{the path where to store the generated plot (\emph{pdf file}). 
           if not specified the plot will be displayed.}
}
\value{
 Returns a list of \code{data.frame} that contain the GO results for each file (or network).
}
\references{
\url{http://david.abcc.ncifcrf.gov/} (DAVID website)

Huang DW, Sherman BT, Lempicki RA. \emph{Systematic and integrative analysis of large gene lists using DAVID Bioinformatics Resources}. Nature Protoc. 2009;4(1):44-57.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{outputGenesPerClusterToDir}}
}
\examples{
    ## get the different datasets path
    petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
    tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  

 \dontrun{

    x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
    ## build the different indexes
    x <- createIndexes(x)
  
    ## build networks connecting each interacting regions
    nets<- buildNetworks(x)

    ## infer the networks
    hlda<- InferNetworks(nets)

    ## Get the list of genes in each cluster by default 
    ## a folder ClustersGenes will be created
    outputGenesPerClusterToDir(hlda,x)

    ## GO enrichment 
    GOEnrich.folder(folder="ClustersGenes/")

}
}

\keyword{methods}
