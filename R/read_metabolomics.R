#' metabolomics read in function
#' This function reads in directory of csv files, using the filenames for method. Can handle multiple methods from multiple patients at the same time
#' Look in given directory for .csv files, read all into a dataframe, adding method name and cleaning up sample ID. Replace NA values with min/10 per compound or with given value
#' Requires dplyr, tidyr, stringr, reshape2
#' @param directory directory containing .csv metabolomic tables in wide format
#' @param na.value if no value present for a compound, replace with this value
#' @param min_over10 when TRUE will override na.value and replace na's with minimum value of the compound / 10
#' @export

read_metabolomics <- function(directory,na.value=1,min_over10=F){
  setwd(dir=directory)
  files <- dir() %>%
    as.data.frame()
  colnames(files) <- "files"

  files <- files %>%
    filter(grepl(".csv",files)) %>%
    mutate(method=gsub("\\.csv","",files),
           method=gsub(" ","\\_",method),
           files=as.character(files))

  data <- list()

  for(i in 1:nrow(files)){
    print(paste0(i," ",files$method[i]))
    data[[i]] <- read.csv(file=files$files[i],header = T,sep=",") %>%
      melt() %>%
      mutate(method=files$method[i])
    colnames(data[[i]]) <- c("compound","sample","peak","method")
  }

  #clean up sample IDs, combine all csv files
  for(i in 1:nrow(files)){
    if(ncol(data[[i]])>4){
      print("ERROR: too many columns... double check the metabolomic data file for accidental entries in blank cells. Unable to combine melted .csv files")
    }else{
    }
  }

  meta <- do.call(rbind,data) %>%
    mutate(Sample_ID=gsub("\\_","",sample),
           Sample_ID=str_extract(pattern="FMT\\.[0-9]{4}[A-Z]+|[0-9]{4}[A-Z]+",sample)) %>%
    #clean up polar methods and metabolite names that start with numeric ID:
    filter(!grepl(compound,pattern="D4")) %>%
    mutate(compound=gsub("[0-9][0-9] ","",compound))

  #find bad sample IDs
  bad_ids <- meta %>%
    filter(is.na(Sample_ID)) %>%
    select(sample) %>%
    unique()

  if(nrow(bad_ids)>0){
    print("Warning: the following Sample_IDs were not formatted correctly, please correct in raw tables or downstream: ")
    print(paste(bad_ids$sample))
  }else{
    print("All sample IDs appear correctly formatted")
  }


  if(min_over10==T){
    #get minimum value for each compound / 10:
    minmethod <- meta %>%
      filter(!is.na(peak),
             !method=="Medications") %>%
      group_by(method) %>%
      summarize(min_method=min(peak)/10)

    #replace na values with min/10:
    meta <- meta %>%
      left_join(minmethod) %>%
      mutate(peak=if_else(is.na(peak),min_method,peak))
  }else{
    meta <- meta %>%
      replace_na(list(peak=na.value))
  }
  return(meta)
}

