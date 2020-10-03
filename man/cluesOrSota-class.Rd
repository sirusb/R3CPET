\name{cluesOrSota-class}
\Rdversion{1.1}
\docType{class}
\alias{cluesOrSota-class}

\title{Wrapper for \code{clues} and \code{sota} S3 classes}
\description{
 This is an S4 virtual union class that defines a new object that can be a \code{\link{sota}}
 or a \code{\link{clues}} class. Now that the \code{clues} method is deprecated, on the \code{\link{sota}} method is supported.
}
\section{Definition}{
\code{setClassUnion("cluesOrSota", c("sota","NULL"))}
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
No methods defined with class "cluesOrSota" in the signature.
}
\references{
Herrero, J., Valencia, A, and Dopazo, J. (2005).\emph{A hierarchical unsupervised growing neural network for clustering gene expression patterns}. Bioinformatics, 17, 126-136.

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
