#' plot metabolomics heatmap
#'
#' @details designed for use with gg.stack with other patient metadata, 16S, etc using 'Sample_ID' on the x axis
#' @param data dataframe containing all metabolite data in long form for one patient. should have the following variables: Sample_ID, compound, peak, method
#' @param row_mean_centered whether to center by mean peak area per compound, default = TRUE
#' @param organize_by either facet by method or order by hierarchical cluster. options are: "method" or "cluster"
#' @param samp_text_size size of Sample_ID label on x axis
#' @examples meta_plot <- plot_meta_hm(meta,row_mean_centered=T,organize_by="cluster")
#' gg.stack(tax,meta.plot,heights=c(1,5))
#' @export

plot_meta_hm <- function(data,row_mean_centered=T,organize_by="method",samp_text_size=18){
  if(row_mean_centered==T){
    #log 2 peak values
    data <- data %>%
      mutate(peak=log(peak,base=2))

    #difference of peak and average, makes column called 'fc' though this isn't exactly fold change
    data <- data %>%
      group_by(compound) %>%
      summarize(avg=mean(peak)) %>%
      left_join(data) %>%
      mutate(fc=peak-avg)
    print("using row mean centered (log2 all values... peak - avg per compound)...")
  }else{
    data <- data %>%
      mutate(fc=peak)
    print("using raw peak areas to draw heatmap...")
  }
  if(organize_by=="method"){
  #metabolomic heatmap (faceted by method, no clustering)
  meta_gg <- data %>%
    ggplot() +
    geom_tile(aes(x=Sample_ID,y=compound,fill=fc)) +
    scale_fill_gradient2(mid="gray92",high="#CD2626",low="navy",midpoint=0) +
    facet_grid(method ~ .,space="free",scales="free_y") +
    theme_bw() +
    theme(strip.text.y=element_text(angle=0,size=14),
          axis.text.x=element_text(size=sampsize,angle=90),
          axis.title.y=element_text(size=18,angle=0)) +
    xlab("Sample_ID") +
    ylab("")
  return(meta_gg)

  }else if(organize_by=="cluster"){

    #cluster code (still very ugly....)
    xx <- data %>%
      dcast(compound ~ Sample_ID,value.var="fc",fun.aggregate=function(x){x[1]})
    row.names(xx) <- xx$compound
    xx <- xx %>%
      select(-compound)
    hc <- hclust(dist(xx,method="euclidean"))
    xx1 <- data.frame(hc$labels,c(1:nrow(as.data.frame(hc$labels))))
    colnames(xx1) <- c("compound","id")
    xx2 <- data.frame(hc$order,c(1:nrow(as.data.frame(hc$order))))
    colnames(xx2) <- c("id","order")
    #finally get order dataframe:
    clust_order <- xx1 %>%
      left_join(xx2)
    #remove intermediates
    rm(xx1,xx2)

    meta_gg <- data %>%
      left_join(clust_order) %>%
      ggplot() +
      geom_tile(aes(x=Sample_ID,y=reorder(compound,order),fill=fc)) +
      scale_fill_gradient2(mid="gray92",high="#CD2626",low="navy",midpoint=0) +
      theme_bw() +
      theme(strip.text.y=element_text(angle=0,size=14),
            axis.text.x=element_text(size=sampsize,angle=90),
            axis.title.y=element_text(size=18,angle=0)) +
      xlab("Sample_ID") +
      ylab("")

    return(meta_gg)

  }else{
    print("Please choose 'method' or 'cluster' for organize_by option...")
  }


}
