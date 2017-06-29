#' Extract Legend Function
#' 
#' This function extracts a legend from a ggplot object, returning a ggplot object


extract_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
