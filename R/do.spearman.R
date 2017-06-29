#' Do Spearman Function
#' 
#' This function lets you perform spearman correlations in dplyr syntax.
#' @param x Numeric Predictor
#' @param y Numeric Outcome
#' @keywords spearman
#' @examples 
#' mtcars %>%
#' group_by(cyl) %>%
#'  summarize(pvalue=as.numeric(do.spearman(wt,mpg)[3]),
#'            rho=as.numeric(do.spearman(wt,mpg)[4]))


do.spearman <- function(x,y) {
  results <- cor.test(x,y,method="spearman")
  return(results)
}