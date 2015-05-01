alpha.Matrix <- function(e) {
  
  prs <- peers(e)
  lmbd <- as.data.frame(lambda(e))
  
  numberOfDMUs <- dim(prs)[1]
  alphaMatrix <- matrix(data = 0, nrow = numberOfDMUs, ncol = numberOfDMUs)
  
  for (i in 1:dim(prs)[2]) {
    alphaMatrix[numberOfDMUs, prs[numberOfDMUs, i]] <- lmbd[numberOfDMUs, i]
  }
  return (alphaMatrix)
}