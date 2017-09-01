#' Rarefy Even Depth 2
#' This function fixes an issue in phyloseq's rarefy even depth function that failed to resample when samples were below seq threshold
#' @export
rarefy_even_depth2 <- function (physeq, sample.size = min(sample_sums(physeq)), rngseed = FALSE,
                                replace = TRUE, trimOTUs = TRUE, verbose = TRUE) {
  rarefaction_subsample <- function (x, sample.size, replace = FALSE) {
    rarvec <- numeric(length(x))
    if (sum(x) <= 0) {
      return(rarvec)
    }
    if (replace) {
      suppressWarnings(subsample <- sample(1:length(x), sample.size, replace = TRUE, prob = x))
    }
    else {
      obsvec <- apply(data.frame(OTUi = 1:length(x), times = x), 1, function(x) {
        rep_len(x["OTUi"], x["times"])
      })
      obsvec <- unlist(obsvec, use.names = FALSE)
      suppressWarnings(subsample <- sample(obsvec, sample.size, replace = FALSE))
    }
    sstab <- table(subsample)
    rarvec[as(names(sstab), "integer")] <- sstab
    return(rarvec)
  }
  if (as(rngseed, "logical")) {
    set.seed(rngseed)
    if (verbose) {
      message("`set.seed(", rngseed, ")` was used to initialize repeatable random subsampling.")
      message("Please record this for your records so others can reproduce.")
      message("Try `set.seed(", rngseed, "); .Random.seed` for the full vector",
              sep = "")
      message("...")
    }
  }
  else if (verbose) {
    message("You set `rngseed` to FALSE. Make sure you've set & recorded\n",
            " the random seed of your session for reproducibility.\n",
            "See `?set.seed`\n")
    message("...")
  }
  if (length(sample.size) > 1) {
    warning("`sample.size` had more than one value. ", "Using only the first. \n ... \n")
    sample.size <- sample.size[1]
  }
  if (sample.size <= 0) {
    stop("sample.size less than or equal to zero. ", "Need positive sample size to work.")
  }
  if (min(sample_sums(physeq)) < sample.size & !replace) {
    rmsamples = sample_names(physeq)[sample_sums(physeq) < sample.size]
    if (verbose) {
      message(length(rmsamples), " samples removed", "because they contained fewer reads than `sample.size`.")
      message("Up to first five removed samples are: \n")
      message(rmsamples[1:min(5, length(rmsamples))], sep = "\t")
      message("...")
    }
    physeq = prune_samples(setdiff(sample_names(physeq),
                                   rmsamples), physeq)
  }
  newsub <- physeq
  if (!taxa_are_rows(newsub)) {
    newsub <- t(newsub)
  }
  newotu <- apply(otu_table(newsub), 2, rarefaction_subsample,
                  sample.size = sample.size, replace = replace)
  rownames(newotu) <- taxa_names(physeq)
  otu_table(newsub) <- otu_table(newotu, TRUE)
  if (trimOTUs) {
    rmtaxa = taxa_names(newsub)[taxa_sums(newsub) <= 0]
    if (length(rmtaxa) > 0) {
      if (verbose) {
        message(length(rmtaxa), "OTUs were removed because they are no longer \n",
                "present in any sample after random subsampling\n")
        message("...")
      }
      newsub = prune_taxa(setdiff(taxa_names(newsub), rmtaxa),
                          newsub)
    }
  }
  if (!taxa_are_rows(physeq)) {
    newsub <- t(newsub)
  }
  return(newsub)
}
