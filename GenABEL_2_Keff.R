keffConv <- function(gwaa.data){
  
  map <- data.frame(SNP.Name = snpnames(gtdata(gwaa.data)))
  map$Order <- 1:nrow(map)
  
  snplist <- data.frame(SNP.Name = snpnames(gtdata(gwaa.data)))
  
  x <- data.frame(All1  = rep(NA, nrow(snplist)), All2  = rep(NA, nrow(snplist)), 
                  Geno1 = rep(NA, nrow(snplist)), Geno2 = rep(NA, nrow(snplist)),
                  Geno3 = rep(NA, nrow(snplist)))
  
  snplist <- cbind(snplist, x)
  
  genotab <- data.frame(t(as.character.gwaa.data(gwaa.data)), stringsAsFactors=F)
  
  genotab <- cbind(data.frame(SNP.Name = row.names(genotab)), genotab)
  
  for(i in 2:ncol(genotab)) genotab[,i] <- gsub("/", "", genotab[,i])
  
  genotab <- merge(snplist, genotab)
  
  
  
  testfunc1 <- function(genrow){
    x <- sort(unlist(unique(c(genrow))))
    y <- substr(x[1], 1, 1)
    y
  }
  
  testfunc2 <- function(genrow){
    x <- sort(unlist(unique(c(genrow))))
    y <- substr(x[length(x)], 2, 2)
    y
  }
  
  genotab$All1 <- apply(genotab[,7:ncol(genotab)], 1, testfunc1)
  genotab$All2 <- apply(genotab[,7:ncol(genotab)], 1, testfunc2)
  
  genotab$Geno1 <- paste(genotab$All1, genotab$All1, sep = "")
  genotab$Geno2 <- paste(genotab$All1, genotab$All2, sep = "")
  genotab$Geno3 <- paste(genotab$All2, genotab$All2, sep = "")
  
  
  for(i in 7:ncol(genotab)){
    
    if(i %in% seq(7,ncol(genotab),10)) print(paste("Individual", i-6))
    genotab[is.na(genotab[, i]), i] <- 0
    genotab[genotab[, i] == genotab[,4], i] <- 1
    genotab[genotab[, i] == genotab[,5], i] <- 2
    genotab[genotab[, i] == genotab[,6], i] <- 3
    genotab[genotab[, i] %in% c("GA", "TC", "CA", "GC", "TA", "TG"), i] <- 2
    
    genotab[, i] <- as.numeric(genotab[, i])
    
  }
  
  head(genotab[,1:10])
  
  genotab <- subset(genotab, select = -c(All1, All2, Geno1, Geno2, Geno3))
  
  genotab <- merge(map, genotab)
  
  genotab <-  genotab[with(genotab, order(Order)), ]
  
  genotab <- subset(genotab, select = -Order)
  
  genotab
  
}