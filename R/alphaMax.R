alpha.Max <- function(alphaMatrix) {
  dimension <- dim(alphaMatrix)[1]
  alphaMax <- matrix(data = 0, nrow = dimension, ncol = 1)
  # Find the max lambda for each DMU
  for (i in 1:dimension){
    alphaMax[i] <- max(alphaMatrix[i,])
  }
  return (alphaMax)
}