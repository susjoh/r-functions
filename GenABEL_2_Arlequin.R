# FOR WINDOWS SYSTEMS ONLY!

genabel2arlequin <- function(gwaa.data, GroupVec, outfile, Title, GroupTogether = T){
  
  genotab <- data.frame(t(as.character.gwaa.data(gwaa.data)), stringsAsFactors = F)
  
  genotab <- t(genotab)
  
  genotab <- data.frame(genotab, stringsAsFactors=F)
  
  pheno <- phdata(gwaa.data)
  
  trait_info <- na.exclude(pheno[,c("id", GroupVec)])
  
  arltab <- NULL
  
  for(g in sort(unique(trait_info[,GroupVec]))){
    
    tempgen <- trait_info[trait_info[,GroupVec] == g,]
    
    space <- matrix(rep("", times = (nsnps(gwaa.data) + 1)),
                    nrow = 1,
                    ncol = (nsnps(gwaa.data) + 1))
    
    arltab <- rbind(arltab, 
                    x <- rbind(cbind(paste("SampleName=\"", GroupVec, "_", g,"\"",sep = ""), space),
                               cbind(paste("SampleSize=", nrow(tempgen), sep = ""), space),
                               cbind("SampleData={", space)))
    
    for(h in tempgen$id){
      
      print(h)
      
      all1 <- NULL
      all2 <- NULL
      
      for(i in 1:ncol(genotab)){
        
        all1 <- cbind(all1, substr(genotab[h,i],1,1)[[1]])
        all2 <- cbind(all2, substr(genotab[h,i],3,3)[[1]])
        
      }
      
      arltab <- rbind(arltab,
                      cbind(h,1,all1),
                      cbind("","",all2))
    }
    
    arltab <- rbind(arltab, cbind("}", space))
  }
  
  for(i in 1:ncol(arltab)) arltab[is.na(arltab[,i]),i] <- "?"
  
  table(arltab[,4])
  
  write.table(arltab, "arltab.temp", quote = F, row.names = F, col.names = F)
  
  outfile <- paste("\"",outfile,"\"", sep = "")
  
  x <- sort(unique(trait_info[,GroupVec]))
  
  system("cmd", input = paste("echo [Profile] > ",outfile, sep = ""))
  system("cmd", input = paste("echo Title=\"", Title,"\" >> ",outfile, sep = ""))
  system("cmd", input = paste("echo NbSamples=", length(x), " >> ",outfile, sep = ""))
  system("cmd", input = paste("echo DataType=STANDARD >> ",outfile))
  system("cmd", input = paste("echo GenotypicData=1 >> ",outfile))
  system("cmd", input = paste("echo LocusSeparator=WHITESPACE >> ",outfile))
  system("cmd", input = paste("echo GameticPhase=0 >> ",outfile))
  system("cmd", input = paste("echo MissingData=\"?\" >> ",outfile))
  system("cmd", input = paste("echo.>> ", outfile))
  system("cmd", input = paste("echo [Data] >> ", outfile))
  system("cmd", input = paste("echo [[Samples]] >> ", outfile))
  
  system("cmd", input = paste("type arltab.temp >> ", outfile))
  
  system("cmd", input = paste("echo.>> ", outfile))
  
  system("cmd", input = paste("echo [[Structure]] >> ", outfile))
  system("cmd", input = paste("echo StructureName=\"A first example of a genetic structure\" >> ", outfile))
  


  
  if(GroupTogether == T){
    system("cmd", input = paste("echo NbGroups = 1 >> ", outfile))
    system("cmd", input = paste("echo Group ={ >> ", outfile))
    
    for(i in 1:length(x)){
      system("cmd", input = paste("echo \"", GroupVec, "_", x[i], "\" >> ", outfile, sep = ""))
    }
    
    system("cmd", input = paste("echo } >> ", outfile))
    system("cmd", input = "rm arltab.temp")
  }
  
  if(GroupTogether == F){
    
    system("cmd", input = paste("echo NbGroups =",length(x)," >> ", outfile))
    
    for(i in 1:length(x)){
      system("cmd", input = paste("echo Group ={ >> ", outfile))
      system("cmd", input = paste("echo \"", GroupVec, "_", x[i], "\" >> ", outfile, sep = ""))
      system("cmd", input = paste("echo } >> ", outfile))
    }
    
    
    system("cmd", input = "rm arltab.temp")
  }
  
}








