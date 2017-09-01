#' old ggstack function
#' ggstack function that still works on older R version
#' @export
gg.stack.el <- function(...,heights = NULL,gg.extras=NULL,gap=0,margin=1,units="inches",as.list=FALSE) {
  grobs <- list(...)
  length.grobs <- length(grobs)
  if (length.grobs<=1) {
    stop("YTError: should have at least 2 grobs")
  }
  top.theme <- theme(plot.margin=unit(c(margin, margin, gap, margin),units),
                     axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  middle.theme <- theme(plot.margin=unit(c(gap, margin, gap, margin),units),
                        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  bottom.theme <- theme(plot.margin=unit(c(gap, margin, margin, margin),units))
  g.top <- grobs[[1]] + top.theme + gg.extras
  g.middle.list <- lapply(grobs[c(-1,-length.grobs)],function(g) {
    g + middle.theme + gg.extras
  })
  g.bottom <- grobs[[length.grobs]] + bottom.theme + gg.extras
  grobs1 <- c(list(g.top),g.middle.list,list(g.bottom))
  grobs2 <- lapply(grobs1,function(g) {
    #gr <- ggplotGrob(g1)
    gr <- ggplotGrob(g)
  })
  nwidths <- max(sapply(grobs2,function(g) length(g$width)))
  grobs3 <- lapply(grobs2,function(g) {
    if (length(g$widths)<nwidths) {
      g <- gtable_add_cols(g,unit(1,"null"))
    }
    return(g)
  })
  max.widths <- do.call(unit.pmax,lapply(grobs3,function(x) x$width))
  grobs4 <- lapply(grobs3,function(g) {
    g$widths <- max.widths
    return(g)
  })
  if (as.list) {
    return(grobs4)
  }
  args <- c(grobs4,list(ncol=1,heights=heights))
  do.call(grid.arrange,args)
}

