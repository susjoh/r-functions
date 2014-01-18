# Functions for populations

PopVecCreator <- function(inputfile) {
  
  system("cmd", input = paste("awk <", inputfile, "\"{print $1}\" > temp.txt"))
  
  system("cmd", input = "sed -n /Pop/{n;p;} < temp.txt > temp2.txt")
  
  x <- read.table("temp2.txt", stringsAsFactors = F)
  
  system("cmd", input = paste("type temp2.txt > ", inputfile, ".poplist", sep = ""))
  
  system("cmd", input = "del temp.txt && del temp2.txt")
  
  return(x)
}


DistmatPipeline <- function(inpvec) {

    system("cmd", input = paste("(for %i in (",paste(inpvec, sep=" ", collapse=" ") ,") do @echo %~i) > bat.txt"))
    
    system("cmd", input = "populations < bat.txt")
    
    system("cmd", input = "del bat.txt")
    
    x <- read.table(inpvec[4], header = F, skip = 3)
    
    return(x)
}


PhyloPipeline <- function(inpvec) {
  
  system("cmd", input = paste("(for %i in (",paste(inpvec, sep=" ", collapse=" ") ,") do @echo %~i) > bat.txt"))
  
  system("cmd", input = "populations < bat.txt")
  
  system("cmd", input = "del bat.txt")
  
  x <-   read.tree(inpvec[length(inpvec)-1], skip = 1)
  
  return(x)
}

    
runPopulations <- function(inputfile1, prefix1, ConductNewAnalysis){
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #   1. Conduct Analyses & read in files   #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #~~ All populations
  
  if(ConductNewAnalysis == TRUE) {
    
    #~~ Determine the population order
    
    poporder1 <- PopVecCreator(inputfile1)
    
    #~~ Distance matrix pipeline:
    input <- c(2, inputfile1, 1, paste(prefix1, ".Dm.dist.txt", sep = ""), 1, "Y")
    distmat1 <- DistmatPipeline(input)
    
    input <- c(2, inputfile1, 1, paste(prefix1, ".Dc.dist.txt", sep = ""), 3, "Y")
    distmat2 <- DistmatPipeline(input)
    
    #~~ Phylogenetic trees
    
    if(length(poporder1) > 2){
      
      # Nei min genetic distance, NJ, 1000 bootstraps
      input <- c(2, inputfile1, 3, 1, 2, 1000, paste(prefix1, "Dm.nj.1000.tre", sep = ""), "Y")
      tree1_wild <- PhyloPipeline(input)
      
      # Nei min genetic distance, NJ
      input <- c(2, inputfile1, 2, 1, 2, paste(prefix1, "Dm.nj.tre", sep = ""), "Y")
      tree2_wild <- PhyloPipeline(input)
      
      # Nei min genetic distance, UPGMA
      input <- c(2, inputfile1, 3, 1, 1, 1000, paste(prefix1, "Dm.upgma.1000.tre", sep = ""), "Y")
      tree3_wild <- PhyloPipeline(input)
      
      
      # Cavalli Sforza, NJ, 1000 bootstraps
      input <- c(2, inputfile1, 3, 3, 2, 1000, paste(prefix1, "Dc.nj.1000.tre", sep = ""), "Y")
      tree4_wild <- PhyloPipeline(input)
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #   2. If analysis has already been done...   #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(ConductNewAnalysis == FALSE) {
    
    poporder1 <- read.table(paste(inputfile1,".poplist", sep = ""), stringsAsFactors = F)
    
    distmat1 <- read.table(paste(prefix1, ".Dm.dist.txt", sep = ""), header = F, skip = 3)
    distmat2 <- read.table(paste(prefix1, ".Dc.dist.txt", sep = ""), header = F, skip = 3)
    
    if(nrow(poporder1) > 2){
      tree1_wild <- read.tree(paste(prefix1, "Dm.nj.1000.tre", sep = ""), skip = 1)
      tree2_wild <- read.tree(paste(prefix1, "Dm.nj.tre", sep = ""), skip = 1)
      tree3_wild <- read.tree(paste(prefix1, "Dm.upgma.1000.tre", sep = ""), skip = 1)
      tree4_wild <- read.tree(paste(prefix1, "Dc.nj.1000.tre", sep = ""), skip = 1)
    }
    
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #   3. Integrate Information from Populations   #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #~~ Sort out the distance matrix
  
  distmat1$V1 <- poporder1$V1
  names(distmat1) <- c("Pop", poporder1$V1)
  
  distmat2$V1 <- poporder1$V1
  names(distmat2) <- c("Pop", distmat2$V1)
  
  #~~ treefiles
  if(nrow(poporder1) > 2){
    newvec <- NULL
    for(i in tree1_wild$tip.label) newvec <- c(newvec, poporder1[poporder1$V2 == i,1])
    tree1_wild$tip.label <- newvec
    
    newvec <- NULL
    for(i in tree2_wild$tip.label) newvec <- c(newvec, poporder1[poporder1$V2 == i,1])
    tree2_wild$tip.label <- newvec
    
    newvec <- NULL
    for(i in tree3_wild$tip.label) newvec <- c(newvec, poporder1[poporder1$V2 == i,1])
    tree3_wild$tip.label <- newvec
    
    newvec <- NULL
    for(i in tree4_wild$tip.label) newvec <- c(newvec, poporder1[poporder1$V2 == i,1])
    tree4_wild$tip.label <- newvec
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #   4. Visual Representations of the Data       #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # trees
  if(nrow(poporder1) > 2){
    plot.phylo(tree1_wild, show.node.label = T, font = 1, cex = 0.8, use.edge.length = F)
    plot.phylo(tree1_wild, show.node.label = T, font = 1, cex = 0.8)
    
    plot.phylo(tree2_wild, show.node.label = T, font = 1, cex = 0.8, use.edge.length = F)
    plot.phylo(tree2_wild, show.node.label = T, font = 1, cex = 0.8)
    
    plot.phylo(tree3_wild, show.node.label = T, font = 1, cex = 0.8, use.edge.length = F)
    plot.phylo(tree3_wild, show.node.label = T, font = 1, cex = 0.8)
    
    plot.phylo(tree4_wild, show.node.label = T, font = 1, cex = 0.8, use.edge.length = F)
    plot.phylo(tree4_wild, show.node.label = T, font = 1, cex = 0.8)
    plot.phylo(tree4_wild, show.node.label = T, font = 1, cex = 0.8, type = "unrooted")
  }
  
  distmat1
  distmat2
}



# 
#     
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # SUPERCEDED



#~~~~~~~~~~~~~~~~~~~~~~~~~
# Populations can be run from the command line as so:
#
#      populations name_of_input_file -"arguments"
#
# Available arguments:
# -phylogeny ind or pop           (default) for phylogenetic trees based on individuals or populations
# -dist method                    (default: Nei standard, Ds) you can choose among: 
#                                      DAS, Dm, Ds, Dc, Da, dmu2, Fst, Cp, Dr, ASD, Dsw, Dr, Dru, Drw, Drl.
# -construct method               (default: upgma) possibilities upgma or nj (Neighbor Joining)
# -bootstrap_ind number           indicate the number of bootstraps to perform on individuals
# -bootstrap_locus number         indicate the number of bootstraps to perform on loci
# -output name_of_treeview_file   indicate the name of the tree file (phylip tree format)
# -level number                   structured populations allows to choose the structuration factor 
#                                     (in the example: town level is 1, building level is 2...).
# 
# example:
#   populations toutc2.txt -phylogeny pop -dist Dm -bootstrap_locus 10000 -output toutc2_10000_Dm.tre

# 
# TreeBuild <- function(inputfile, prefix, distance, construct, poporder, bootstrap_locus = 0) {
#   
#   require(ape)
#   
#   outfile <- paste(prefix,".", distance, ".", construct, bootstrap_locus, ".tre", sep = "")
#   
#   if(bootstrap_locus == 0) {
#     
#     cmd1 <- paste("populations ", inputfile, 
#                   " -phylogeny pop -dist ", distance,
#                   " -construct ", construct,
#                   " -output ", outfile, sep = "")
#   }
#   
#   if(bootstrap_locus != 0) {
#     
#     cmd1 <- paste("populations ", inputfile, 
#                   " -phylogeny pop -dist ", distance,
#                   " -construct ", construct,
#                   " -bootstrap_locus ", bootstrap_locus,
#                   " -output ", outfile, sep = "")
#   }
#   
#   system("cmd", input = cmd1)
#   
#   newpops <- paste("POP_", 1:length(poporder), ":", sep = "")
#   
#   system("cmd", input = "mkdir tmp")
#   
#   for(i in 1:length(newpops)){
#     system("cmd", input = paste("sed -e \"s/",newpops[i],"/", poporder[i],":/g\" < ", outfile, " > tmp/", outfile, sep = ""))
#     system("cmd", input = paste("mv tmp/", outfile, " .", sep = ""))
#   }
#   
#   system("cmd", input = "rmdir tmp")
#   
#   treex <- read.tree(outfile, skip = 1)
#   
#   return(treex)
#   
# }