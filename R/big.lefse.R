#' Big Lefse function
#'
#' Convert melted taxonomy file into Lefse formatted .txt file
#'
#' @export
big.lefse <- function(t,dictionary,filename="galaxy_lefse.txt"){
  xx <- t %>% group_by(sample,Kingdom) %>% summarise(pctseqs=sum(pctseqs))
  colnames(xx) <- c("sample","taxonall","pctseqs")
  xx <- as.data.frame(xx)
  xx2 <- t %>% group_by(sample,Kingdom,Phylum) %>% summarise(pctseqs=sum(pctseqs))
  xx2$taxonall <- paste(xx2$Kingdom,xx2$Phylum,sep="|")
  xx2 <- as.data.frame(xx2) %>% select(sample,taxonall,pctseqs)
  xx3 <- t %>% group_by(sample,Kingdom,Phylum,Class) %>% summarise(pctseqs=sum(pctseqs))
  xx3$taxonall <- paste(xx3$Kingdom,xx3$Phylum,xx3$Class,sep="|")
  xx3 <- as.data.frame(xx3) %>% select(sample,taxonall,pctseqs)
  xx4 <- t %>% group_by(sample,Kingdom,Phylum,Class,Order) %>% summarise(pctseqs=sum(pctseqs))
  xx4$taxonall <- paste(xx4$Kingdom,xx4$Phylum,xx4$Class,xx4$Order,sep="|")
  xx4 <- as.data.frame(xx4) %>% select(sample,taxonall,pctseqs)
  xx5 <- t %>% group_by(sample,Kingdom,Phylum,Class,Order,Family) %>% summarise(pctseqs=sum(pctseqs))
  xx5$taxonall <- paste(xx5$Kingdom,xx5$Phylum,xx5$Class,xx5$Order,xx5$Family,sep="|")
  xx5 <- as.data.frame(xx5) %>% select(sample,taxonall,pctseqs)
  xx6 <- t %>% group_by(sample,Kingdom,Phylum,Class,Order,Family,Genus) %>% summarise(pctseqs=sum(pctseqs))
  xx6$taxonall <- paste(xx6$Kingdom,xx6$Phylum,xx6$Class,xx6$Order,xx6$Family,xx6$Genus,sep="|")
  xx6 <- as.data.frame(xx6) %>% select(sample,taxonall,pctseqs)
  xx7 <- t %>% group_by(sample,Kingdom,Phylum,Class,Order,Family,Genus,Species)
  xx7$taxonall <- paste(xx7$Kingdom,xx7$Phylum,xx7$Class,xx7$Order,xx7$Family,xx7$Genus,xx7$Species,sep="|")
  xx7 <- as.data.frame(xx7) %>% select(sample,taxonall,pctseqs)
  t1 <- rbind(xx,xx2,xx3,xx4,xx5,xx6,xx7)
  sss <- dcast(data=t1,formula=sample~taxonall,value.var="pctseqs",fill=0,fun.aggregate=sum)
  ssss <- left_join(dictionary,sss)
  lefse <- ssss
  lefse1 <- as.data.frame(t(lefse))
  lefse2 <- t(lefse)
  #lefse2 <- lefse2
  lefse2 <- lefse2[-1,]
  write.table(lefse2, filename, sep = "\t",col.names=FALSE,quote=FALSE)
  getwd()
}
