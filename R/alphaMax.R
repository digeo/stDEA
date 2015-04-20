alpha.Max <- function(alphaMatrix) {
  alphaMax <- matrix(data = 0, nrow = dim(alphaMatrix)[1], ncol = 1)
  # Find the max lambda for each DMU
  for (i in minimumDMUs:numberOfDMUs){
    alphaMax[i] <- max(alphaMatrix[i,])
  }
  return alphaMax
}