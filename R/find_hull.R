#' find hull function
#'
#' this function finds convex hull to draw polygons onto 2d pcoa plots per group
#' @examples
#' #first create ordination object (ord) using phyloseq
#' df <- as.data.frame(as.matrix(ord$vectors))
#' df$sample <- row.names(df)
#' df <- df %>%
#'   left_join(samp)
#'
#'   hulls <- ddply(df,"area",fund_hull)
#'
#'   #sample plot:
#'   df %>%
#' ggplot(aes(x=Axis.1,y=Axis.2)) +
#'   geom_point(size=4,aes(shape=day,color=group)) +
#'   theme_bw(base_size=18) +
#'   coord_equal() +
#'   geom_polygon(data=hulls,aes(fill=area),alpha=.2) +
#'   scale_color_manual(values=c("MNV"="steelblue2","Strep"="salmon")) +
#'   scale_fill_manual(values=c("MNV"="steelblue2","Strep"="salmon","1"="chartreuse4","-3"="gray33")) +
#'   xlab("Axis.1 [27.7%]") +
#'   ylab("Axis.2 [20.1%]") +
#'   ggtitle("PCoA Bray")
#' @export
find_hull <- function(df) {df[grDevices::chull(df$Axis.1,df$Axis.2), ]}
