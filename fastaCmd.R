fastaCmd <- function(db, seqIDs, newfasta){
  
  newfile <- "temp0.txt"
  write.table(data.frame(droplevels(unique(seqIDs))), newfile, quote = F, row.names = F, col.names = F)
  
  system("cmd", input = paste("fgrep -f", newfile,"--line-number", db, "> temp1.out"))
  system("cmd", input = paste("grep \">\" --line-number", db, "> temp2.out"))
  
  seqline <- read.table("temp1.out", sep = ":")
  for(i in 1:length(seqIDs)){
    if(i == 1) system("cmd", input = paste("sed -n \"/", seqIDs[i],"/{n;p;}\" < temp2.out > temp3.out", sep = ""))
    if(i != 1) system("cmd", input = paste("sed -n \"/", seqIDs[i],"/{n;p;}\" < temp2.out >> temp3.out", sep = ""))
  }
  
  refline <- read.table("temp3.out", sep = ":")
  
  combline <- rbind(data.frame(seqline = seqline$V1, ref = "start"), 
                    data.frame(seqline = refline$V1, ref = "stop"))
  
  combline <-  unique(combline[with(combline, order(seqline)), ])
  
  for(i in seq(1,nrow(combline),2)){
    if(i == 1) system("cmd", input = paste("sed -n \"", combline$seqline[i], ",", (combline$seqline[i + 1]-1), "p\" ", db, " > ", newfasta, sep = ""))
    if(i != 1) system("cmd", input = paste("sed -n \"", combline$seqline[i], ",", (combline$seqline[i + 1]-1), "p\" ", db, " >> ", newfasta, sep = ""))
  }

  system("cmd", input = "del temp0.txt temp1.out temp2.out temp3.out")
  
  }