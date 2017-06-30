#' Humann2 Output Read In
#' 
#' Read in all output files for any number of samples into a list of dataframes, containing gene abundances, pathway abundances, and pathway coverage for each sample
#' 
#' Output information:
#' data[[1]] = genefam (gene family abundance table)
#' data[[2]] = pathabu (pathway abundance table)
#' data[[3]] = pathcov (pathway coverage table)
#' 
#' Requires pbapply and stringr libraries
#' @examples 
#' a <- "/Volumes/castoricenter/Eric.Littmann/U01/Processed_Data/Humann2/"
#' data <- humann2_readin(a)
#' 


humann2_readin <- function(file_dir){
  #requires plyr and dplyr, load plyr first
  setwd(file_dir)
  print("reading in humann2 files...")
  humann2process <- function(input) {
    input$tax <- str_extract(input[,1],"\\|.+")
    input$tax <- gsub("\\|","",input$tax)
    tax <- input$tax
    xx <- as.data.frame(tax) %>% separate(tax,into=c("Genus","Species"),"\\.",fill="left")
    yy <- cbind(input,xx)
    yy$tax[is.na(yy$tax)] <- "Total Community"
    yy$sample <- colnames(yy)[2]
    yy$type <- colnames(yy)[1]
    colnames(yy)[2] <- "Abundance"
    colnames(yy)[1] <- "Variable"
    return(yy)
  }
  
  #coerces all data into a big list of dataframes:
  
  data <- ldply(pblapply(pblapply(dir(),read.delim),humann2process),data.frame) %>%
    mutate(type=gsub("X..","",type))
  print("creating list of dataframes...")
  #return(data)
  datalist <- NULL
  datalist[[1]] <- data %>% filter(type=="Gene.Family")
  path <- data %>% filter(type=="Pathway") %>% 
    mutate(subtype=str_extract(sample,"Abundance|Coverage"),
           pathwayname=str_extract(Variable,".+\\:"),
           pathwayname=gsub("\\:","",pathwayname),
           category=str_extract(Variable,
                                "biosynthesis and salvage|biosynthesis|degradation|fermentation|salvage|
                                cycle|resistance|oxidation|reduction|metabolism|recycling|editing|interconversion|
                                processing|superpathway"))
  path$category[is.na(path$category)] <- "other"
  datalist[[3]] <- path %>% filter(subtype=="Coverage")
  datalist[[2]] <- path %>% filter(subtype=="Abundance")
  print("[[1]] = genefam (gene family abundance table)")
  print("[[2]] = pathabu (pathway abundance table)")
  print("[[3]] = pathcov (pathway coverage table)")
  return(datalist)
}
