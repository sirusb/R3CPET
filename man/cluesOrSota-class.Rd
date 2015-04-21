\name{cluesOrSota-class}
\Rdversion{1.1}
\docType{class}
\alias{cluesOrSota-class}

\title{Wrapper for \code{clues} and \code{sota} S3 classes}
\description{
 This is an S4 virtual union class that defines a new object that can be a \code{\link{sota}}
 or a \code{\link{clues}} class.
}
\section{Definition}{
\code{setClassUnion("cluesOrSota", c("clues", "sota","NULL"))}
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
No methods defined with class "cluesOrSota" in the signature.
}
\references{
Herrero, J., Valencia, A, and Dopazo, J. (2005).\emph{A hierarchical unsupervised growing neural network for clustering gene expression patterns}. Bioinformatics, 17, 126-136.
Wang, S., Qiu, W., and Zamar, R. H. (2007). \emph{CLUES: A non-parametric clustering method based on local shrinking. Computational Statistics & Data Analysis}, Vol. 52, issue 1, pages 286-298.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}

\seealso{
\code{\link{ChromMaintainers}}
}
\examples{
showClass("cluesOrSota")
}
\keyword{classes}
