\name{loadPPI-methods}
\docType{methods}

\alias{loadPPI-methods}
\alias{loadPPI,ChiapetExperimentData-method}
\alias{loadPPI}

\title{ Setting the background protein-protein interaction}
\description{
 This method enables the user to define the PPI network as a background network. 
 The user can provide his own PPI or use the \code{HPRD} or the \code{Biogrid} PPI
 incorporated in the package.
}
\usage{
\S4method{loadPPI}{ChiapetExperimentData}(object,type=c("HPRD","Biogid"),customPPI= NULL, 
              filter=FALSE, term ="GO:0005634", annot=NULL, 
              RPKM= NULL, threshold=1 )
}
\arguments{
\item{object}{
  a \code{\link{ChiapetExperimentData}} to load the PPI into as an \code{\link{igraph}} object
}
\item{type}{
 if \code{customPPI = NULL}, this parameter indicats which of the available PPI to use \code{HPRD} or \code{Biogrid} (\emph{build 3.2.100}). 
 by default, the \code{HPRD} network is used.
}
\item{customPPI}{
 If the user wants to use his own PPI interaction network (for example for another species), he can provide an 
 \code{\link{igraph}} object or a path to an \code{"ncol"} formatted graph (a file with two columns indicating the interacting parts
 and an optional third column to indicate the weight). It is preferred that the user provides an \code{igraph} object, to avoid
 any problems when parsing the \code{"ncol"} file.
}
\item{filter}{
 This parameter indicates whether the user want to filter the selected PPI or not. \code{filter = FALSE} means that the PPI (provided by the user or by the \code{HPRD} and \code{Biogrid} networks) will be used as is. 
 if \code{filter = TRUE},the proteins in the PPI will be filtered according to their position in cell ( by default proteins located 
 in the nucleus are kept and the other removed).
 
 In addition to the filtering by location, if the user wants just to keep the proteins that show a certain amount of expression 
 he can provide a gene expression table to the \code{RPKM} parameter and set the \code{threshold} for filtering.
}
\item{term}{
 The GO term used to for filtering. By default, only protein that are located in the nucleus are kept (\code{term = "GO:0005634"}).
 If the user want to use another annotation table he can pass a \code{\link{data.frame}} object to the parameter \code{annot}.
}
\item{annot}{
 If the user wants to provide his own annotation data-set, he can pass a \code{\link{data.frame}} object to this parameter.
 The passed \code{data.frame} should have at least two columns named respectively: \code{geneSymbol} that contains the gene names, 
 and \code{cellular_component_term} that contain the term.
  \preformatted{
    geneSymbol                    cellular_component_term
    FAU Ribosome      (GO:0005840); Nucleolus (GO:0005730)
    ALDH3A1   Cytoplasm (GO:0005737); Nucleus (GO:0005634)
    ASCL1     Nucleus (GO:0005634); Cytoplasm (GO:0005737)
  }
}
\item{RPKM}{
 A \code{\linkS4class{data.frame}} that contains the expression of each gene in the PPI. The \code{\linkS4class{data.frame}} should 
 at least have 2 columns. The first one contains the gene symbol (should be the same as the one used in the PPI) and the second 
 gives the expression.
}
\item{threshold}{
 Threshold value used to filter gene expression. All genes with expression value less than \code{threshold} are removed.
}
}
\value{
A \code{\link{ChiapetExperimentData}} object in which the \code{ppi} slot is populated as an \code{\link{igraph}} object filtered 
according to the specified conditions.
}
\references{
Prasad, T. S. K. et al. (2009) \emph{Human Protein Reference Database - 2009 Update}. Nucleic Acids Research. 37, D767-72.


Chatr-Aryamontri A, Breitkreutz BJ et al. \emph{The BioGRID Interaction Database: 2013 update}. Nucleic Acids Res. 2012 Nov 30


M.N Djekidel et al,\emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, in press, 2015
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{ChiapetExperimentData}}, \code{\link{loadTFBS}} , \code{\link{loadPETs}}, \code{\link{createIndexes}}
}
\examples{

## Create a ChiapetExperimentData object
x <- ChiapetExperimentData(ppiType= "HPRD")

## Loading the default HPRD network without filtering

x <- loadPPI(x,type="HPRD")
ppi(x)

## Using the HPRD network and filtering using the GO:0005634
x <- loadPPI(x,type="HPRD", filter=TRUE)
ppi(x)

data(RPKMS)
x <- loadPPI(x,type="HPRD",filter=TRUE,annot= NULL, RPKM= RPKMS, threshold = 5)
ppi(x)
}


\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
