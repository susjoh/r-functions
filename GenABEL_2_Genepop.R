# gwaa.object <- tenogen_unlinked
# GroupVec <- "SeaAge1v3"
# outfile <- "test.txt"

#~~ Last update 2013-03-08

genabel2genepop <- function(gwaa.object, GroupVec = NA, outfile){
  
  if(is.na(GroupVec)){
    gwaa.object <- add.phdata(gwaa.object, rep("Pop1", nids(gwaa.object)), "TempVec")
    GroupVec <- "TempVec"
  }
  
  x <- as.character.snp.data(gtdata(gwaa.object))
  x <- data.frame(x, stringsAsFactors = F)
  
  x <- cbind(row.names(x), x)
  names(x)[1] <- "id"
  
  x[is.na(x)] <- "0000"
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "A", replacement = "01")); print("Converted A's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "C", replacement = "02")); print("Converted C's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "G", replacement = "03")); print("Converted G's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "T", replacement = "04")); print("Converted T's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "/", replacement = ""))
  
  #head(x[,1:10])
  
  #~~ Create data frame to start
  
  pheno <- phdata(gwaa.object)
  
  x <- merge(x, pheno[,c("id", GroupVec)])
  
  # reorder the columns
  
  genepop <- data.frame(Population = paste(GroupVec, x[,GroupVec], sep = "_"))
  genepop$sep <- rep(",",times = nrow(genepop))
  genepop <- cbind(genepop, x[,2:ncol(x)-1])
  
  #head(genepop[,1:10])
  
  system("cmd",input = paste("echo Title line: GENEPOP ANALYSIS > ", outfile))
  system("cmd",input = paste("echo ID > ", outfile, ".ids", sep = ""))
  
  write.table(names(genepop)[4:ncol(genepop)], outfile, row.names = F, quote = F, col.names = F, append=T)
  
  for(i in sort(unique(genepop$Population))) {
    system("cmd", input = paste("echo Pop >> ", outfile))
    x <- subset(genepop, Population == i)
    write.table(data.frame(ID = x$id), paste(outfile, ".ids", sep = ""), row.names = F, col.names = F, sep = "\t", quote = F, append = T)
    write.table(x[,-3], outfile, row.names = F, col.names = F, sep = "\t", quote = F, append = T)
  }
  
}

