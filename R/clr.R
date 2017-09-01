#' CLR transformation function
#' Performs centered log ratio transformation to normalize count data
#' @examples
#' t <- get.otu.melt(phy)
#' t %>%
#' group_by(sample) %>%
#' summarize(pctseqs=clr(pctseqs))
#' @export
gm_mean <- function(x, na.rm=TRUE){
# The geometric mean, with some error-protection bits
exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

clr <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}
