
require(msm)

estVQ <- function(asremlmodel, dataframe, animalid, snpid){
  
  #~~ estimate allele frequencies
  freqs <- table(unique(dataframe[,c(animalid, snpid)])[,2])
  if(length(freqs) == 3){
    p <- (freqs[1] + 0.5*freqs[2])/sum(freqs)
    q <- 1-p
  }
  if(length(freqs) != 3) stop("not enough genotypes")
  
  #~~ estimate a and d
  
  fixeftab <- summary(asremlmodel, all = T)$coef.fixed
  fixeftab <- fixeftab[grep(snpid, row.names(fixeftab)),]
  
  a = (fixeftab[1,1] - fixeftab[3,1])/2
  # if(a < 0) a <- a * -1
  
  d = a + fixeftab[2,1]
  
  Vq <- 2*p*q*(a + d*(q - p))^2
  Va <- summary(asremlmodel, all = T)$varcomp[paste("ped(", animalid, ", var = T)!ped", sep = ""),]$component
  VarExplained <- Vq/(Vq + Va)
  
  C <- MCMCglmm::Tri2M(asremlmodel$Cfixed, FALSE)
  diag(C)
  sqrt(diag(C))   # match s.e. of models (have a look0
  
  C<-C[2:3,2:3]
  C
  
  x1 <- C[1,1]   # sampling variance of the het effect
  x2 <- C[2,2]   # sampling variance of the hom effect
  beta <- summary(asremlmodel, all = T)$coef.fixed[2:3, 1]
  
  X <- 2*p*q
  Y <- q^2
  
  Vq.se <- deltamethod(~X*(-x2/2 + (-x2/2 + x1)*Y)^2, beta, C)  # standard error
  
  results <- list(Vq, Va, VarExplained, Vq.se)
  names(results) <- c("Vq", "Va", "VarExplained", "Vq.se")
  
  results
  
}
