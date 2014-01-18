# grm.object <- grm.gcta
# ids.object <- ids.gcta 
# phenoframe <- horn.data.oldrm  
# idcol <- 1

makeGRM <- function(grm.object, ids.object, phenoframe, idcol = 1){
  
  names(phenoframe)[idcol] <- "ANIMAL"
  phenoframe$ANIMAL <- as.factor(phenoframe$ANIMAL)
  
  elements<-length(ids.object[,1])
  X <- diag(elements)
  
  X[upper.tri(X, diag=TRUE)] <- grm.object[,4]
  X.grm <- X + t(X) - diag(diag(X))
  nam<-as.factor( ids.object[,2]);
  
  #THE FOLLWING PART IS ESSENTIAL TO MAKE THE MATRIX CONTAIN THE EXACT SAME IDS AS THE PHENOTYPE FILE
  rownames(X.grm)<-nam ;
  colnames(X.grm)<-nam
  
  attr(X.grm,"rowNames")<-as.factor(rownames(X.grm))
  
  
  d<- intersect(ids.object$V2,phenoframe$ANIMAL)#VECTOR OF IDs shared between phenotype and genotype file
  b<-match(phenoframe$ANIMAL,d,nomatch=20000)
  phenoframe<-as.data.frame(subset(phenoframe,b!=20000) )    #PHENOTYPE FILE THAT MACHES THE GENOTYPE FILE
  
  cn<-match(ids.object[,2],d,nomatch=20000)
  nam2<-as.vector(subset(nam,cn!=20000))
  
  rc<-match(rownames(X.grm),d,nomatch=20000)
  X.grm<-X.grm[rc!=20000,]    #DELETE REDUNDANT ROWS
  cc<-match(colnames(X.grm),d,nomatch=20000)
  X.grm<-X.grm[,rc!=20000]   #DELETE REDUNDANT COLUMNS
  X2.grm<-ginv(X.grm)     #CREATE GENERALIZED INVERSE MATRIX
  rownames(X2.grm)<-nam2
  colnames(X2.grm)<-nam2
  attr(X2.grm,"rowNames")<-as.factor(nam2)
  
  
  phenoframe$ANIMAL<-as.factor(phenoframe$ANIMAL)
  
  x <- list()
  x[["grm"]] <- X2.grm
  x[["phenoframe"]] <- phenoframe
  
  x
}





