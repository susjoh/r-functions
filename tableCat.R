tableCat <- function(inFrame) {
    outText <- paste(names(inFrame), collapse = " | ")
    outText <- c(outText, paste(rep("---", ncol(inFrame)), collapse = " | "))
    invisible(apply(inFrame, 1, function(inRow) {
        outText <<- c(outText, paste(inRow, collapse = " | "))
    }))
    return(outText)
}