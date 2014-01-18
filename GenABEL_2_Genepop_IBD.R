Create_Genepop_IBD <- function(gwaa.data, fileprefix){
  
  #~~ create the genepop file
  
  x <- as.character.snp.data(gtdata(gwaa.data))
  x <- data.frame(x, stringsAsFactors = F)
  
  x <- cbind(row.names(x), x)
  names(x)[1] <- "id"
  
  x[is.na(x)] <- "000000"
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "A", replacement = "001")); print("Converted A's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "C", replacement = "002")); print("Converted C's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "G", replacement = "003")); print("Converted G's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "T", replacement = "004")); print("Converted T's")
  x[,2:ncol(x)] <- as.data.frame(sapply(x[,2:ncol(x)], gsub, pattern = "/", replacement = ""))
  
  head(x[,1:10])
  
  #~~ add cumu.dist information
  
  distinfo <- phdata(tenogen_QCpass.1to6)[,c("id", "cumu.dist")]
  
  distinfo$DummyDist <- "0,"
  x1 <- merge(distinfo, x)
  head(x1[,1:10])
  
  
  #~~ sort by cumulative distance
  
  x1 <- x1[with(x1, order(cumu.dist)), ]
  
  
  #~~ create the file
  
  IDlist <- x1$id
  write.table(IDlist, paste(fileprefix, ".IDlist.txt", sep = ""), row.names = F, quote = F, col.names = F)
  
  system("cmd", input = paste("echo Title line: GENEPOP ANALYSIS > ", fileprefix, ".txt", sep = ""), show.output.on.console=F)
  
  write.table(names(x1)[4:ncol(x1)],"tempout2.txt", row.names = F, quote = F, col.names = F)
  system("cmd",input = "type tempout2.txt >> test.txt", show.output.on.console=F)
  
  x1$id <- "Popxxx"
  write.table(x1,"tempout2.txt", row.names = F, quote = F, col.names = F)
  system("cmd", input = "type tempout2.txt >> test.txt", show.output.on.console=F)
  system("cmd", input = "rm tempout2.txt", show.output.on.console=F)
  print("sed starting")
  system("cmd", input = paste("sed -e \"s/xxx /\\n/g\" < test.txt >> ", fileprefix, ".txt", sep = ""), show.output.on.console=F)
  system("cmd", input = "rm test.txt", show.output.on.console=F)
  print("Finished!")
  
}