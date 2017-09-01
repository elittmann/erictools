#'General 16S Abundance Plot
#'
#'This function reads uparse output files from the specified directory and generates a 16S abundance plot at the species level
#'
#' @param dir Directory containing the uparse output files (otu table in biom format, fasta of otu reps, taxonomy table, phylogenetic tree). Should be string ie "~/Downloads"
#' @param threshold Minimum number of sequences required for labeling taxons at species level
#' @param taxtextsize Font size of labeled taxons
#' @param labeltextsize Font size of sample labels on the x axis
#' @export
el_abundance_plot <- function(dir,threshold=100,taxtextsize=3,labeltextsize=12){
  setwd(dir)

  #import uparse data:
  biom.file <- "total.8.otu-tax.biom"
  seq.file <- "total.5.repset.fasta"
  tree.file <- "total.10.tree"
  tax.file <- "total.5.repset.fasta.blastn.refseq_rna.txt"

  biom <- import_biom(biom.file)
  seq <- import_qiime(refseqfilename=seq.file)
  tree <- read.tree.uparse(tree.file)
  tax <- read.blastn.file(tax.file) %>% set.tax()
  phy <- merge_phyloseq(biom,seq,tree)
  tax_table(phy) <- tax
  #remove clutter:
  rm(biom.file,seq.file,tax.file,tree.file,biom,seq,tax,tree)

  samp <- get.samp(phy) %>%
    mutate(pool=str_extract(pattern="\\.\\..+",sample),
           pool=gsub("\\.\\.","",pool))

  pools <- unique(samp$pool)

  #sample_data(phy) <- set.samp(samp)

  t <- get.otu.melt(phy)
  pal <- get.yt.palette(t)

  gg <- t %>%
    arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    mutate(Species = factor(Species, levels = unique(Species))) %>%
    group_by(sample) %>%
    arrange(Species) %>%
    mutate(cum.pct = cumsum(numseqs), y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>%
    ungroup() %>%
    dplyr::select(-cum.pct) %>%
    mutate(tax.label = ifelse(numseqs >= threshold, as.character(Species), ""),
           tax.label=gsub(" ", "\n",tax.label)) %>%
    ggplot() +
    geom_bar(aes(x=sample,y=numseqs,fill=Species),stat="identity",color="black") +
    theme(legend.position="none",
          axis.text.x=element_text(angle=90,size=labeltextsize)) +
    scale_fill_manual(values=pal) +
    geom_text(aes(x=sample,y=y.text,label=tax.label),lineheight=.6,size=taxtextsize) +
    ylab("Number of Sequences") +
    xlab("") +
    ggtitle(paste(pools))

  return(gg)


}
