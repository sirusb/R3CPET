\name{InferNetworks-methods}
\docType{methods}
\alias{InferNetworks}
\alias{InferNetworks-methods}
\alias{InferNetworks,NetworkCollection-method}
\title{ Network construction using Hierarchical Dirichlet Process}
\description{
 This methods applies a Hierarchical Dirichlet Process (HDP) algorithm on the collection of proteins networks to infer 
 the set of chromatin loop-maintainer proteins.
 HDP are non-parametric Bayesian models widely used in document classification as it enables us to model
 datasets with a mixtures of classes.
 In our case,  we suppose that different kinds of networks are involved in maintaining the different loops. 
 Thus, to make an analogy with topic modeling, we each DNA-interaction maintaining protein network as a document
 and each edge in this network as word. Thus, the task is to say which word (\emph{edge}) belongs to which 
 topic (\emph{chromatin-maintainer family}).  
 The method implementation is based on the C++ code of Chong Wang and David Blei with adaptation to Rcpp and removal of the 
 dependency on the Gnu Scientific Library.
 }
\usage{
\S4method{InferNetworks}{NetworkCollection}(object,thr =0.5,max_iter = 500L, max_time = 3600L, \dots)
}
\arguments{
 \item{object}{ 
 a \code{\link{NetworkCollection}} object in which the list of protein interactions associated with each DNA interaction 
 is already populated.
 }
 \item{thr}{ Used to select the top protein interaction in each inferred chromatin-maintainer family.
 In HDP each topic (\emph{Chromatin-maintainer family}) is considered as a distribution over words (\emph{edges}), 
 thus, for each topic we consider the words that capture \code{threshold} percent of the topic to be the top words.
 For example, in topic1, we first rank edges by partnership probability to topic1 in a decreasing order, and we take the top edges 
 that capture 50\% of the partnership. }
 \item{max_iter}{maximum number of iterations (befault 500).}
 \item{max_time}{maximum runing time (3600 sec).}
  \item{\dots}{ You can pass additional paramters to control the behaviour of the HDP model. The possible paramters are \emph{eta}, \emph{alpha} and \emph{gamma}. \emph{eta} controls how edges are assigned to CMNs on the global level. smaller \emph{eta} values will lead to sparce edge-to-CMN assignment, which \emph{eta} >1 leads to more uniform assignments. \emph{gamma} on the other hand, controls the number of CMNs, smaller \emph{gamma} values will produce a small number of CMNs and \emph{gamma}>1 will favor the generation of more. \emph{alpha} controls the sparcity at the local PPI. smaller \emph{alpha} value force edges to be conrolled by a small number of CMNs, while lagrger values leads to more uniform distribution. By default \emph{eta = 0.01}, \emph{gamma =1} and \emph{alpha =1}. 
  }    
}
\value{
 Returns a \code{\link{ChromMaintainers}} object that contains the list of inferred networks and the probability
 of each edge in each network.
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
\code{\link{NetworkCollection}}, \code{\link{ChromMaintainers}} 
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
    hlda
}
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
