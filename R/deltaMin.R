delta.min <- function(deltaMatrix) {
  dimension <- dim(deltaMatrix)[1]
  deltaMin <- matrix(data = 0, nrow = dimension, ncol = 1)
  
  # Find the DeltaMin
  for (i in 1:dimension){
    deltaMin[i] <- 0
    
    for (j in 1:i){
      if((deltaMatrix[i, j] != M) && (deltaMatrix[i, j] > deltaMin[i])){
        deltaMin[i] <- deltaMatrix[i, j]
      }
    }
    if(deltaMin[i] == 0){
      deltaMin[i] <- 1
    }
  }
  return (deltaMin)
}