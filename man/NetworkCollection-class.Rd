\name{NetworkCollection-class}
\docType{class}
\alias{NetworkCollection}
\alias{NetworkCollection-class}

\alias{networks}
\alias{networks,NetworkCollection-method}
\alias{sizes}
\alias{sizes,NetworkCollection-method}
\alias{TF}
\alias{TF,NetworkCollection-method}

\title{protein interaction networks maintaining DNA loops}
\description{
 The class \code{NetworkCollection} stores information about the set of protein networks that maintains DNA interactions.
}
\details{
 The \code{NetworkCollection} contains three main information:
 \emph{(i)} the set of edges in each network maintaining each DNA loop, 
  \emph{(ii)} the number of edges in each network and \emph{(iii)}
  the set of TF involved in all the networks.
 
}
\usage{
 NetworkCollection(networks, sizes,  TFCollection)
}
\arguments{
    \item{networks}{a \code{list} of \code{list}, each list contains the set of edges in each network}
    \item{sizes}{the sizes of each network. The should correspond to the sizes of the networks}
    \item{TFCollection}{the set of all the TF involved in all the interactions.}
}
\section{Accesors}{
  \describe{
    \item{networks}{gets the list of networks}
    \item{sizes}{gets the vector containing the size of each network }
    \item{TF}{gets the list of involved TF (after filtering)}
	 }
}

\value{
  a \code{NetworkCollection} object.
}
\references{
Mohamed Nadhir D, Yang C et al \emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, ....
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}

\seealso{
\code{\link{InferNetworks}}, \code{\link{ChiapetExperimentData}}, \code{\link{buildNetworks}}
}
\examples{
showClass("NetworkCollection")
}
\keyword{classes}
