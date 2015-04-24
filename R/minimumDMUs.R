minimumDMUs <- function(numberOfInputs, numberOfOutputs) {
  return (max(numberOfInputs * numberOfOutputs, 3 * (numberOfInputs + numberOfOutputs)))
}