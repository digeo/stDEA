delta.Matrix <- function(alphaMatrix) {
  
  # Define the M
  M <- 100000
  
  dimension <- dim(alphaMatrix)[1]
  deltaMatrix <- matrix(data = M, nrow = dimension, ncol = dimension)
  
  # Calculating the Delta matrix
  for (i in 1:dimension){
    for (j in 1:dimension){
      if (alphaMatrix[i,j] > 0){
        deltaMatrix[i,j] <- (i - j)
      }
    }
  }
  return (deltaMatrix)
}