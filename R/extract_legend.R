#' Extract Legend Function
#' This function extracts a legend from a ggplot object, returning a ggplot object
#' @export
extract_legend<-function(ggplot_object){
  tmp <- ggplot_gtable(ggplot_build(ggplot_object))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
