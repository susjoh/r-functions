
readGenepop <- function(inputfile){
  
  
  x <- read.table(inputfile, sep="\t", fill=T, skip = length(readLines(inputfile))-1)
  if(x[1,2] == ",")  nloci <- ncol(x) - 2
  if(x[1,2] != ",")  nloci <- ncol(x) - 1
  
  loci <- read.table(inputfile, skip = 1, nrows=nloci, stringsAsFactors = F)
  
  
  x <- read.table(inputfile, sep="\t", fill=T, skip = nloci + 2, stringsAsFactors=F)
  
  head(x)[,1:5]
  head(x)[,(ncol(x) - 4):ncol(x)]
  
  if(x[1,2] == ",") x <- x[,-2]
  
  names(x) <- c("ID", loci$V1)
  
  
  x[,1] <- gsub(",", "", x[,1])
  
  popdivs <- c(0, grep("Pop", x[,1], ignore.case=T))
  
  x$Pop <- NA
  for(i in 1:length(popdivs)){
    start <- popdivs[i]
    if(i != length(popdivs)) end   <- popdivs[i + 1] -1  
    if(i == length(popdivs)) end   <- nrow(x)
    
    x$Pop[start:end] <- paste("Pop", i, sep = "")
  }
  
  x <- x[-grep("Pop", x[,1], ignore.case=T),]

  
  
  return(x)
  
}