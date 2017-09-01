#' otu names function
#' this function adds numbers to the front of otu numbers to keep them unique when combining phyloseq objects from different runs
#' @examples
#' #need to give unique OTU names in otu and tax table before combining phyloseq objects from 2 mothur runs!
#'
#'newotu <- otunames("0000",otu_table(mothurphy2))
#'newtax <- otunames("0000",tax_table(mothurphy2))
#'mothurphy2 <- phyloseq(newotu,newtax)
#'phy <- merge_phyloseq(mothurphy,mothurphy2)
#' @export
otunames <- function(x,n) {
  row.names(x) <- paste0(n,row.names(x))
  return(x)
}
