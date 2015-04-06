\name{EnsemblToHGNC}
\alias{EnsemblToHGNC}
\docType{methods}
\title{
Ensemble to HGNC conversion
}
\description{
This helper method uses the \code{biomaRt} package to convert 
Ensembl ids to HGNC ids.
}
\usage{
EnsemblToHGNC(EnsemblIDs)
}
\arguments{
 \item{EnsemblIDs}{ a \code{character} vector with Ensembl IDs.
 }
}
\value{
returns a \code{\link{data.frame}} containing the Ensembl ID and his corresponding HGNC gene id and Name 
plus a description of the gene.
}
\author{
Mohamed Nadhir Djekidel (\email{nde12@mails.tsinghua.edu.cn})
}
\seealso{
\code{\link{EntrezToHGNC}}
}
\examples{
\dontrun{
EnsemblIDs<-c("ENSG00000164548","ENSG00000118515","ENSG00000105705",
        "ENSG00000177414","ENSG00000108179")

EnsemblToHGNC(EnsemblIDs)
}
}
\keyword{methods}
