\name{ChiapetExperimentData-class}
\docType{class}
\alias{ChiapetExperimentData-class}
\alias{ChiapetExperimentData}

\alias{pet}
\alias{pet-methods}
\alias{pet,ChiapetExperimentData-method}
\alias{pet<-}
\alias{pet<-,ChiapetExperimentData-method}
\alias{tfbs}
\alias{tfbs-methods}
\alias{tfbs,ChiapetExperimentData-method}
\alias{tfbs<-}
\alias{tfbs<-,ChiapetExperimentData-method}
\alias{ppi}
\alias{ppi-methods}
\alias{ppi,ChiapetExperimentData-method}
\alias{ppi<-}
\alias{ppi<-,ChiapetExperimentData-method}

\title{3CPET used raw data}
\description{
 The \code{ChiapetExperimentData} class is a container for storing the set of raw data used by 3CPET to do the prediction. 
}
\details{
 The \code{ChiapetExperimentData} class stores the genomic coordinates of the 
 ChIA-PET interactions, the binding sites of the different transcription factor (TFBS) 
 and the background protein-protein interaction (PPI) network used to infer the final 
 chromatin maintainer networks.
}
\usage{
ChiapetExperimentData(pet='', tfbs='', ppi=NULL,  
                      ## loadPETs options
                      IsBed=TRUE, petHasHeader=FALSE, dist=1000,
                      ## loadTFBS options
                      tfbsHasHeader=FALSE,
                      ## loadPPI options
                      ppiType=c("HPRD","Biogid"),
                      filter=FALSE, term="GO:0005634", annot=NULL,
                      RPKM= NULL, threshold=1
                      )
}

\arguments{
    \item{pet}{(\emph{optional}) a ChIAP-PET interactions file path or a \code{\link{GRanges}} object. 
    The \code{GRanges} object should have a column named \code{PET_ID}.    
    details on the file format can be found on the \code{\link{loadPETs}} help page.}
    
    \item{tfbs}{(\emph{optional}) a file path to a BED file containing transcription factors binding site or a \code{GRanges} object. The \code{GRanges} object should have a metadata column named \code{TF}.}
      
    \item{ppi}{(\emph{optional}) an \code{igraph} object that contains a user defined protein-protein interaction network.
    if this parameter is not specified, the \code{ppiType} paramter will be used.}
    
    \item{IsBed}{(\emph{optional}) considered only if the \code{pet} parameter is a file path. 
    More info about this paramters can be found in the \code{\link{loadPETs}} help page. } 
    
    \item{petHasHeader}{(\emph{optional}) \code{logical}. Indicates if the ChIA-PET interactions file has a header or not.} 
    
    \item{dist}{(\emph{optional}) The size of the region to consider arround the center of the interacting regions.}
    
    \item{tfbsHasHeader}{(\emph{optional}) \code{logical}. Indicates if the TFBS file has a header or not.} 
    
    \item{ppiType}{(\emph{optional}) considered only if the \code{ppi} paramter is not specified. 
    This paramter tell the pakage to load one of the PPI (\code{HPRD} or \code{Biogrid}) shiped with the package. }
    
    \item{filter}{(\emph{optional}) \code{logical}. whether of not to filter the PPI network.
      if the \code{RPKM} paramter is specified then the RPKM dataset incorporated with the package will be used.
      if you want to to your own way of filtering, you ca set \code{filter = FALSE} and pass an already processed PPI
      to the \code{ppi} paramterer.}
    
    \item{term}{(\emph{optional}) the GO term used to filter the nodes of the PPI. this is different from the \code{filter} parameter.
    in the \code{filter} parameter the PPI nodes are filtered by their gene expression, while in the \code{term} parameter
    they are filtered by their genomic location. by default \code{"GO:0005634"} is used for filtering.}
    
    \item{annot}{ (\emph{optional}) the annotation dataset used for filtering by default the \code{\link{geneLocations}} is used.
    The user can pass a custom \code{data.frame}. For more details check the \code{\link{loadPPI}} help page.}
    
    \item{RPKM}{(\emph{optional}) a \code{data.frame} object that contains the expression value of each genes.
    by default the \code{\link{RPKMS}} dataset will be used (expression value in K562 celline).
    For more information about the format of the data passed to this parameter please check the \code{\link{loadPPI}}}
    \item{threshold}{ (\emph{optional}) threshold value used to filter gene expression. default: 1. }
      
}

\section{Slots}{
  \describe{
    \item{pet}{:
    Object of class \code{\linkS4class{GRanges}} that stores the genomic coordinated of the interactions.
    it can be populated using the method \code{\link{loadPETs}}}
    \item{tfbs}{:
    Object of class \code{\linkS4class{GRanges}} that stores the TF binding site.
    it can be populated using the method \code{\link{loadTFBS}}.
    NOTE: the TFBS locations can be obtained from a ChIP-Seq experiment or a motif finding software.
    for more information on the format of the provided data check \code{\link{loadTFBS}}}
    \item{ppi}{
    Object of class \code{"igraph"} used as the background PPI for further analysis.
    it can be populated using the method \code{\link{loadPPI}} }
    \item{.dt}{
    Object of class \code{"list"} contains a collection of \code{data.table} 
    serving as indexes used internally by the package (not expected to be manipulated by the user).
    it can be populated using the method \code{\link{createIndexes}} }
  }
}

\section{Accessors}{
 The following methods can be used to get the content of a \code{ChiapetExperimentData} object \code{x} :
 \describe{
    \item{}{
      \code{pet(x), pet(x) <- value}:
      Get ChIA-PET interactions encoded as a \code{\linkS4class{GRanges}} object in \code{x}. The returned \code{\linkS4class{GRanges}} objects 
      contains an attribute \code{PET_ID} in which the left side have an id of the form \code{PET#\\d+\\.1}
      and the right side interaction have an id of the form \code{PET#\\d+\\.2}. for more information check \code{\link{loadPETs}}
      \preformatted{
          seqnames             ranges strand |      PET_ID
             <Rle>          <IRanges>  <Rle> | <character>
      [1]     chr1 [1240734, 1242734]      * |     PET#1.1
      [2]     chr1 [1242224, 1244224]      * |     PET#1.2
      [3]     chr1 [1282208, 1284208]      * |     PET#2.1
      [4]     chr1 [1283334, 1285334]      * |     PET#2.2
      [5]     chr1 [1370371, 1372371]      * |     PET#3.1
      [6]     chr1 [1371822, 1373822]      * |     PET#3.2       
  }
    }
    \item{}{
      \code{tfbs(x), tfbs(x) <- value}:
      Get the \code{\linkS4class{GRanges}} storing the transcription factor binding sites.
    }
    \item{}{
      \code{ppi(x), ppi(x) <- value}:
      Returns an \code{\link{igraph}} object used as a background PPI. check the \code{\link{loadPPI}} for more information.
    }
  }
}

\value{
Constructs a \code{\link{ChiapetExperimentData}} object with the specified fields populated.
}

\references{
Li G, Fullwood MJ, Xu H et al.\emph{ChIA-PET tool for comprehensive chromatin interaction analysis with paired-end tag sequencing}. Genome Biology 2010, 11(2):R22

Mohamed Nadhir D, Yang C et al \emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, ....
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}


\seealso{
\code{\link{loadPETs}}, \code{\link{loadTFBS}} , \code{\link{loadPPI}}
}
\examples{
## for example Reading ChIA-PET interaction results generated from ChIA-PET tool
## it should be formatted as follow:

## -------------------------------------------------------------------------------------------
## chromleft startleft endleft chromright startright endright counts      pvalue      qvalue
##      chr1    872113  879175       chr1     933836   938416     12 1.84529e-30 6.90983e-28
##      chr1    874165  879175       chr1     933340   938306     10 1.23139e-25 3.58932e-23
##      chr1    889676  896594       chr1     933897   938982     13 4.91311e-36 2.33753e-33
##      chr1    898753  907581       chr1     931133   939571     19 0.00000e+00 0.00000e+00
##      chr1    910103  918775       chr1     930834   938627     15 2.20004e-43 1.32812e-40
##      chr1    919314  922154       chr1     934212   937864      6 3.70292e-21 7.88551e-19  
##---------------------------------------------------------------------------------------------

## The counts, pvalue and qvalue fields are not considered in our case 
## it is up to the user to filter the interactions.

## The TFBS should be a BED file that contain the chromosome, start, end and the TF name 

\dontrun{

  ## load the different datasets
  petFile <- file.path(system.file("example",package="R3CPET"),"HepG2_interactions.txt")  
  tfbsFile <- file.path(system.file("example",package="R3CPET"),"HepG2_TF.txt.gz")  
  
  x <- ChiapetExperimentData(pet = petFile, tfbs=  tfbsFile, IsBed = FALSE, ppiType="HPRD", filter= TRUE) 
  ## build the diffrent indexes
  x <- createIndexes(x)
  x
  
  ## Pass objects instead of files.
  pet <- sample(pet(x),size = 20,replace = TRUE )
  tfbs <- sample(tfbs(x), size=300, replace=TRUE)
  ppi <- ppi(x)
  tst <- ChiapetExperimentData(pet = pet, tfbs= tfbs, ppi=ppi) 
  tst <- createIndexes(tst)
  tst
}
}
\keyword{classes}
