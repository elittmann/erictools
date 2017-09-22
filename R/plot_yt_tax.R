#'General 16S Abundance Plot
#' @details This function reads uparse output files from the specified directory and generates a 16S abundance plot at the species level
#' @param phylo phyloseq object (can use read_uparse to get uparse outputs into phyloseq...)
#' @param pctseqs When TRUE will return plots with % abundance, when false will return plots using number of sequences not normalized to 1
#' @param threshold Minimum number of sequences (numseqs) or % abundance (pctseqs) required for labeling taxons at species level
#' @param taxtextsize Font size of labeled taxons
#' @param labeltextsize Font size of sample labels on the x axis
#' @export
plot_yt_tax <- function(phylo,pctseqs=T,threshold=.1,taxtextsize=3,labeltextsize=12){

  t <- get.otu.melt(phy)
  pal <- get.yt.palette(t)

  if(pctseqs==F){
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
    xlab("")

  return(gg)
  }else{
    gg <- t %>%
      arrange(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
      mutate(Species = factor(Species, levels = unique(Species))) %>%
      group_by(sample) %>%
      arrange(Species) %>%
      mutate(cum.pct = cumsum(pctseqs), y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>%
      ungroup() %>%
      dplyr::select(-cum.pct) %>%
      mutate(tax.label = ifelse(pctseqs >= threshold, as.character(Species), ""),
             tax.label=gsub(" ", "\n",tax.label)) %>%
      ggplot() +
      geom_bar(aes(x=sample,y=pctseqs,fill=Species),stat="identity",color="black") +
      theme(legend.position="none",
            axis.text.x=element_text(angle=90,size=labeltextsize)) +
      scale_fill_manual(values=pal) +
      geom_text(aes(x=sample,y=y.text,label=tax.label),lineheight=.6,size=taxtextsize) +
      ylab("Number of Sequences") +
      xlab("")

    return(gg)

 }
}
