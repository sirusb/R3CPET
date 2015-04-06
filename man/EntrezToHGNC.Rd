\name{EntrezToHGNC}
\alias{EntrezToHGNC}
\docType{methods}
\title{
Entrez to HGNC conversion
}
\description{
This helper method uses the \code{biomaRt} package to convert 
Entrez ids to HGNC icS.
}
\usage{
EntrezToHGNC(EntrezID)
}
\arguments{
 \item{EntrezID}{ a \code{character} vector with Entrez IDs.
 }
}
\value{
returns a \code{\link{data.frame}} containing the Entrez ID and his corresponding HGNC gene id and Name 
plus a description of the gene.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{EnsemblToHGNC}}
}
\examples{
\dontrun{
EntrezID <-c("2114","9757","5886","9373","6921",
        "4088","7006","6196","10054","10945")

EntrezToHGNC(EntrezID)
}
}
\keyword{methods}