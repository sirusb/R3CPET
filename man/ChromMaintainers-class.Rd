\name{ChromMaintainers-class}
\Rdversion{1.1}
\docType{class}
\alias{ChromMaintainers}
\alias{ChromMaintainers-class}

%accesors
%\alias{annotateExpression,ChromMaintainers,data.frame-method}
%\alias{cluster,ChromMaintainers-method}
\alias{getClusters}
\alias{getClusters-methods}
\alias{getClusters,ChromMaintainers-method}

% methods
\alias{GenerateGmlNetworks}
\alias{GenerateGmlNetworks-methods}
\alias{GenerateGmlNetworks,ChromMaintainers-method}
%\alias{outputGenesPerClusterToDir,ChromMaintainers,ChiapetExperimentData-method}
%\alias{getRegionsIncluster,ChromMaintainers,ChiapetExperimentData,numeric-method}
\alias{networks,ChromMaintainers-method}
%\alias{plot3CPETRes,ChromMaintainers-method}
\alias{topEdges}
\alias{topEdges-methods}
\alias{topEdges,ChromMaintainers-method}
\alias{topNodes}
\alias{topNodes-methods}
\alias{topNodes,ChromMaintainers-method}
%\alias{visualizeCircos,ChromMaintainers,ChiapetExperimentData,numeric-method}

\title{Chomatin maintainer networks}
\description{
 The  \code{ChromMaintainers} holds information about the inferred network by the method \code{\link{InferNetworks}}. 
 It contains the list of inferred networks as \code{\link{igraph}} object, a list of edges and a list of proteins.
 In addition to an \code{\link{HLDAResult}} object that contains the final probabilities calculated by the HLDA algorithm.
}
\usage{
 ChromMaintainers( maintainers,topEdges,topNodes, clusRes = NULL, networks = list())
}
\arguments{
    \item{maintainers}{Object of class \code{"HLDAResult"} that contains the HLDA results. }
    \item{topEdges}{a \code{"matrix"} containing the list the top edges per network. }
    \item{topNodes}{a \code{"matrix"} containing the list the top proteins per network.}
    \item{networks}{the list of infered networks as an \code{\link{igraph}} objects list}
    \item{clusRes}{Object of class \code{"cluesOrSota"} describing the assignment of each DNA-interaction 
    to a chromatin-maintainer network according to their enrichment.}
}
\section{Accesors}{
  if \code{x} is a \code{ChromMaintainers} object the following accessors can be applied :
  \describe{
    \item{\code{networks(x)}}{gets the list of networks as \code{\link{igraph}} objects}
    \item{\code{topNodes(x)}}{gets a \code{matrix} object that contains the list of top proteins per network }
    \item{\code{topEdges(x)}}{gets a \code{matrix} object that contains the list of top proteins per network}
    \item{\code{getClusters(x)}}{returns the clustering results of DNA-interaction into groups according to their partnership 
    enrichment profile to the set of inferred chromatin maintainer networks.}
	 }
}

\section{Methods}{
 Many plotting and annotation methods are associated with this class.
  \describe{
    \item{\code{annotateExpression(object, RPKMS)}}{To add the gene expression attribute to the \code{\link{igraph}} objects}
    \item{\code{clusterInteractions(object, method, nbClus )}}{
     To cluster the DNA-interactions according to their partnership enrichment profile.}
    \item{\code{GenerateGmlNetworks(object,...)}}{Creates the list of \code{\link{igraph}} object from the \code{topEdges} matrix}
    \item{\code{outputGenesPerClusterToDir(hdaRes,data,path, ...)}}{ 
    get the list genes belonging to each DNA-interaction cluster.}
    \item{\code{getRegionsIncluster(hdaRes,data, cluster, ...)}}{returns the coordinates of the DNA interactions for a given cluster.}
    \item{\code{GOEnrich.networks(object, pval=0.05, GOlimit= 5,path="")}}{
    do a \code{GO} enrichment of the elements of each inferred network.}
    \item{\code{plot3CPETRes(object, path, W, H, type,byEdge, netPerRow, layoutfct, ...)}}{
         provide different type of plots to visualize the results}
    \item{\code{visualizeCircos(object, data, cluster, chrLenghts)}}{ Draws a circos plot of the DNA interactions in a given cluster.}
	 }
}

\value{
	a \code{ChromMaintainers} object.
}

\references{
\url{https://www.cs.princeton.edu/~blei/topicmodeling.html} (C. Wang's hdp code)

Chong Wang, John Paisley and David M. Blei, \emph{Online variational inference for the hierarchical Dirichlet process}
.In AISTATS 2011

Mohamed Nadhir D, Yang C et al \emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, ....
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{InferNetworks}}, \code{\link{ChromMaintainers}} , \code{\link{HLDAResult}}
}
\examples{
showClass("ChromMaintainers")
}
\keyword{classes}
