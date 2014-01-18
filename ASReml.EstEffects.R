ASReml.EstEffects <- function(model){
  x <- summary(model)$varcomp
  totvar   <- sum(x[,1])
  x$Effect <- x$gamma/totvar
  x$SE <- NA
  
  object <- model
  pframe <- as.list(object$gammas) 
  
  denominator <- "1 "
  
  for(effname in names(pframe)[1:length(pframe)-1]){
    denominator <- paste(denominator, "+ `", effname, "` ", sep = "")
  }
  
  denominator <- paste("(", denominator, ")", sep = "")
  
  for(effname in names(pframe)[1:length(pframe)-1]){
    transform <- eval(parse(text = paste("`placeholder` ~ `", effname, "`/", denominator, sep = "")))
    
    tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe) 
    X <- as.vector(attr(tvalue, "gradient")) 
    tname <- if (length(transform) == 3) {
      transform[[2]] 
    } else "" 
    V <- object$ai 
    n <- length(pframe) 
    i <- rep(1:n, 1:n) 
    if (length(V) != length(i)) 
      stop("vcov matrix incompatible with\nnumber of variance components") 
    j <- sequence(1:n) 
    k <- 1 + (i > j) 
    se <- sqrt(sum(V * X[i] * X[j] * k))
    x[effname,"SE"] <- se
  }
  
  x
}