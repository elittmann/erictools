#' read_uparse
#'
#' uparse readin function
#'
#' read outputs from Ying's uparse pipeline into a phyloseq object
#' @param directory directory containing uparse output files (usually called 'uparse')
#' requires phyloseq, ape, yingtools2, dplyr
#' @examples phy <- read_uparse("Z:/Eric/uparse")
#' @export

read_uparse <- function(directory){
  setwd(directory)
  biom.file <- "total.8.otu-tax.biom"
  seq.file <- "total.5.repset.fasta"
  tree.file <- "total.10.tree"
  tax.file <- "total.5.repset.fasta.blastn.refseq_rna.txt"


  biom <- import_biom(biom.file)
  seq <- import_qiime(refseqfilename=seq.file)
  tree <- read.tree.uparse(tree.file)
  taxa_names(tree) <- gsub("\\_\\_","\\;",taxa_names(tree))

  tax <- read.blastn.file(tax.file) %>% set.tax()
  phy <- merge_phyloseq(biom,seq,tree)
  tax_table(phy) <- tax

  rm(biom.file,seq.file,tax.file,tree.file,biom,seq,tax,tree)
  return(phy)

}
