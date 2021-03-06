% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rarefy_even_depth2.R
\name{rarefy_even_depth2}
\alias{rarefy_even_depth2}
\title{Rarefy Even Depth 2
This function fixes an issue in phyloseq's rarefy even depth function that failed to resample when samples were below seq threshold}
\usage{
rarefy_even_depth2(physeq, sample.size = min(sample_sums(physeq)),
  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
}
\description{
Rarefy Even Depth 2
This function fixes an issue in phyloseq's rarefy even depth function that failed to resample when samples were below seq threshold
}
