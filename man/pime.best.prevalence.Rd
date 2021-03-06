% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pime.best.prevalence.R
\name{pime.best.prevalence}
\alias{pime.best.prevalence}
\title{Best Prevalence}
\usage{
pime.best.prevalence(prev.list, variable, method.dist = "bray", ...)
}
\arguments{
\item{prev.list}{List phyloseq objects with the calculated prevalences for each interval. The output of
pime.prevalence()}

\item{variable}{Any variable in the metadata to be analyzed}

\item{method.dist}{Distance/dissimilarity. Default is Bray-Curtis. Can be other dissimilarity distance supported by vegan::vegdist().}

\item{...}{Aditional parameters passed to ranger::ranger() on random forest classification.}
}
\description{
This function is the core of PIME. It runs a Permanova analysis on distance/dissimilarity matrix,
and performs a classification of variables with random forest, for each prevalence interval returned by
pime.prevalence(). The output is a table with the Permanova results, OOB error rate and
the number of remaining OTUs and sequences after prevalence filtering.
}
\examples{
phylist=pime.split.by.variable(restroom, "Environment")
prev=pime.prevalence(phylist)
pime.best.prevalence(prev, "Environment", method.dist="bray")
}
\keyword{OOB}
\keyword{permanova}
\keyword{prevalence}
