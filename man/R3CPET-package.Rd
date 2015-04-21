\name{R3CPET-package}
\alias{R3CPET-package}
\alias{R3CPET}
\docType{package}
\title{
3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process
}
\description{
The main goal of 3CPET is to try to infer the set of protein networks that are likely to  work together inorder to maintain chromatin loops obtained by a ChIA-PET experiment.
It is based on an idea silimar to the one used for document classification. It starts first by building a PPI network for each chromatin interaction, then uses an HDLA (Hierarchical Dirichlet Latent Allocation) model to infer the set of networks that are enriched together.

}
\details{
\tabular{ll}{
Package: \tab R3CPET\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-23\cr
License: \tab GPL (>= 3.0)\cr
}
}
\author{
Written by M.N.Djekidel
Maintainer: Mohamed Nadhir Djekidel <nde12@mails.tsinghua.edu.cn>
}
\references{
M.N Djekidel et al,\emph{3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process}, in press, 2015
}
\keyword{ package }
\seealso{
\code{\link{ChiapetExperimentData}}, \code{\link{ChromMaintainers}} , \code{\link{HLDAResult}}
}
