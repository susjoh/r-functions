ASReml.LRT <- function(model1, model2) {
  
  # check that DF are the same for each model
  if(summary(model1)$nedf != summary(model2)$nedf) stop("d.f. mismatch: Models are not directly comparable")
  
  chi.stat <- 2*diff(sort(c(summary(model1)$loglik, summary(model2)$loglik)))
  
  pchisq(chi.stat, 1, lower.tail=F)
  
}
