#'Best Prevalence
#'
#'This function is the core of PIME. It runs a Permanova analysis on distance/dissimilarity matrix,
#'and performs a classification of variables with random forest, for each prevalence interval returned by
#'pime.prevalence(). The output is a table with the Permanova results, OOB error rate and
#'the number of remaining OTUs and sequences after prevalence filtering.
#'
#'
#'@param prev.list List phyloseq objects with the calculated prevalences for each interval. The output of
#'pime.prevalence()
#'@param variable Any variable in the metadata to be analyzed
#'@param method.dist Distance/dissimilarity. Default is Bray-Curtis. Can be other dissimilarity distance supported by vegan::vegdist().
#'@param ... Aditional parameters passed to ranger::ranger() on random forest classification.
#'@keywords prevalence permanova OOB
#'@examples phylist=pime.split.by.variable(restroom, "Environment")
#' prev=pime.prevalence(phylist)
#' pime.best.prevalence(prev, "Environment", method.dist="bray")
#' @importFrom phyloseq "otu_table"
#' @importFrom phyloseq "sample_data"
#'@export
pime.best.prevalence = function (prev.list, variable, method.dist="bray",...) {
  set.seed(2125)
  perm=list()
  randon=list()
  gs <- as(object = sample_data(prev.list[[1]]), Class = "data.frame")
  Variable = as.factor(gs[, variable])
  for (i in prev.list){
    if ((phyloseq::taxa_are_rows(i)==TRUE)==TRUE){
      train=t(otu_table(i))
    } else {
    train=otu_table(i)}
    d=vegan::vegdist(train, method = method.dist)
    adonis=vegan::adonis(d ~ Variable, data = gs)
    perm[[length(perm)+1]]=adonis[[1]]
    print("Calculating...")
    response = Variable
    # Combine them into 1 data frame
    training.set <- data.frame(response, train)
    train.model = ranger::ranger(response ~ ., data = training.set,...)
    randon[[length(randon)+1]]=train.model$prediction.error
    print("Done")
  }
  #names tables from Lista as the names of the tables inside list.core
  names(perm) <- paste("Prevalence", names(prev.list))
  #gets only the first line, all columns of every table inside Lista
  results1 <- do.call(rbind, lapply(perm,`[`,1,))
  Interval= paste(as.numeric(names(prev.list)), "%", sep = "")
  results2 <- do.call(rbind, randon)
  Nseqs=sapply(prev.list, function(z) sum(phyloseq::sample_sums(z)))
  OTUs=sapply(prev.list, phyloseq::ntaxa)
  OOB.err=cbind(results1,results2,Interval,OTUs,Nseqs)
  return(OOB.err)
}
