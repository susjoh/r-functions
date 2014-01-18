

genabel2plink <- function(gwaa.data, phenoname, recode12 = TRUE, outname){
  
  require(car)
  
  genotab <- data.frame(t(as.character.gwaa.data(gwaa.data)), stringsAsFactors = F)
  
  genotab <- t(genotab)
  
  genotab <- data.frame(genotab, stringsAsFactors=F)
  
  genotab2 <- matrix(nrow = nrow(genotab), ncol = ncol(genotab)*2)
  
  for(i in 1:ncol(genotab)){
    genotab2[,((2*i)-1)] <- substr(genotab[,i], 1, 1)
    genotab2[,(2*i)]     <- substr(genotab[,i], 3, 3)
    
    if(i %in% seq(0, 100, ncol(genotab))) print(paste(i, "Loci Converted"))
  }
  
  genotab2[is.na(genotab2)] <- 0
  
  genotab2 <- data.frame(genotab2)
  
  genotab2$id <- row.names(genotab)
  
  pheno <- phdata(gwaa.data)
  pheno$sex[pheno$sex == 0] <- 2
  
  if(recode12 == TRUE) pheno[,phenoname] <- recode(pheno[,phenoname], '0 = 1; 1 = 2; NA = 0')
  
  pedfile <- data.frame(Family = pheno$id,
                        ID     = pheno$id,
                        Father = rep(0, times = length(pheno$id)),
                        Mother = rep(0, times = length(pheno$id)),
                        Sex    = pheno$sex,
                        Pheno  = pheno[,phenoname])
  
  pedfile$ID == genotab2$id
  
  pedfile <- cbind(pedfile, genotab2)
  
  pedfile <- pedfile[,-ncol(pedfile)]
  
  
  write.table(pedfile, paste(outname, ".ped", sep = ""), row.names = F, col.names = F, quote = F)
  
  mapfile <- data.frame(SNP = names(genotab),
                        Order = 1:length(names(genotab)))
  
  mapfile <- merge(mapfile, data.frame(Chromosome = chromosome(gwaa.data),
                                       SNP = snp.names(gwaa.data),
                                       cM = rep(0, times = length(snp.names(gwaa.data))),
                                       Position = map(gwaa.data)), all.x = T)
  
  mapfile <-   mapfile[with(mapfile, order(Order)), ]
  mapfile <- subset(mapfile, select = -Order)
  mapfile <- mapfile[,c(2,1,3,4)]
  
  mapfile[is.na(mapfile)] <- 0

  write.table(mapfile, paste(outname, ".map", sep = ""), row.names = F, col.names = F, quote = F)
  
}







