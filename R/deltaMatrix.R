delta.Matrix <- function(alphaMatrix) {
  
  # Define the M
  M <- 100000
  
  dimetion <- dim(alphaMatrix)[1]
  deltaMatrix <- matrix(data = M, nrow = dimetion, ncol = dimetion)
  
  # Calculating the Delta matrix
  for (i in 1:dimetion){
    for (j in 1:dimetion){
      if (alphaMatrix[i,j] > 0){
        deltaMatrix[i,j] <- (i - j)
      }
    }
  }
  return (deltaMatrix)
}