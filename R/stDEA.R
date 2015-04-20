#' Spatio-Temporal DEA efficiency
#' 
#' Estimates a ST-DEA frontier and calculates efficiency measures a la Farrell.
#' @param X Inputs of firms to be evaluated, a K x m matrix of observations of K firms with m inputs (firm x input). In case TRANSPOSE=TRUE the input matrix is transposed to input x firm.
#' @param Y Outputs of firms to be evaluated, a K x n matrix of observations of K firms with n outputs (firm x input). In case TRANSPOSE=TRUE the output matrix is transposed to output x firm.
#' @param RTS Text string or a number defining the underlying DEA technology / returns to scale assumption.
#' @param ORIENTATION Input efficiency "in" (1), output efficiency "out" (2), and graph efficiency "graph" (3). For use with DIRECT , an additional option is "in-out" (0).
#' @return The results are returned in a Farrell object with the following components. The last three components in the list are only part of the object when SLACK=TRUE.
#' 
#' @export
stDEA <- function(X, Y, RTS = "vrs", ORIENTATION = "out"){
  
  rts <- c("fdh", "vrs", "drs", "crs", "irs", "irs2", "add", "fdh+", "fdh++", "fdh0")
  
  if (missing(RTS)) 
    RTS <- "vrs"
  if (is.numeric(RTS)) {
    RTStemp <- rts[1 + RTS]
    RTS <- RTStemp
  }
  RTS <- tolower(RTS)
  if (!(RTS %in% rts)) 
    stop(paste("Unknown scale of returns:", RTS))
  
  orientation <- c("in", "out")
  
  if (is.numeric(ORIENTATION)) {
    ORIENTATION_ <- orientation[ORIENTATION + 1]
    ORIENTATION <- ORIENTATION_
  }
  
  ORIENTATION <- tolower(ORIENTATION)
  
  if (!(ORIENTATION %in% orientation)) {
    stop(paste("Unknown value for ORIENTATION:", ORIENTATION))
  }
  
  # Define the M
  M <- 100000
  
  # Step
  stp <- 0.01
  
  repeats <- 1 / stp
  
  numberOfDMUs <- dim(X)[1]
  numberOfOutputs <- dim(Y)[2]
  numberOfInputs <- dim(X)[2]
  
  minimum.DMUs <- max(numberOfInputs * numberOfOutputs, 3 * (numberOfInputs + numberOfOutputs))
  stop.DMU <- minimum.DMUs
  
  if(minimumDMUs > numberOfDMUs){
    stop("Not enough Decision Making Units, you should have at least : ", minimumDMUs)
  }
  
  efficiency.dea <- rep(NA, times = minimumDMUs - 1)
  
  alphaMatrix <- matrix(data = 0, nrow = numberOfDMUs, ncol = numberOfDMUs)
  
  while(stop.DMU <= number.Of.DMUs){
    tempX <- X[1:stopDMU,]
    tempY <- Y[1:stopDMU,]
    e <- dea(tempX, tempY, RTS, ORIENTATION)
    efficiency.dea <- c(efficiency.dea, eff(e)[stopDMU])
    prs <- dea.peers(e)
    lmbd <- as.data.frame(lambda(e))
    for (i in 1:dim(prs)[2]) {
      alphaMatrix[stopDMU, prs[stop.DMU, i]] <- lmbd[stop.DMU, i]
    }
    stop.DMU <- stop.DMU + 1
  }
  
  # Define the constraintTypes
  constraintTypes <- matrix(data = ">=", nrow = (numberOfOutputs + 1), ncol = 1)
  
  # Define the rhs
  rightHandSide <- matrix(data = 0, nrow = numberOfOutputs, ncol = 1)
  
  # ST-DEA
  for (i in minimumDMUs:numberOfDMUs){
    objectiveFunction <- matrix(data = 0, nrow = i, ncol = 1)
    
    
    for(Wsp in seq(0, 1, stp)){
      Wt <- 1 - Wsp
      
      lpmodel <- make.lp((numberOfOutputs + 2), (i + 1))
      
      # The 1st column, is the column for the phi
      set.column(lpmodel, 1, c((-Y[i,]), 1, 0))
      
      for (j in 1:i){
        objectiveFunction[j] <- (((Wsp / alphaMax[i]) * alphaMatrix[i, j]) - ((Wt / deltaMin[i]) * deltaMatrix[i, j]))
        set.column(lpmodel, (j + 1), c((Y[j,]), 0, 1))
      }
      
      set.objfn(lpmodel, c(0, objectiveFunction[1:i]))
      set.constr.type(lpmodel, c(constraintTypes[1:(numberOfOutputs + 1)], "="))
      set.rhs(lpmodel, c(rightHandSide[1:numberOfOutputs], 1, 1))
      
      set.type(lpmodel, 2:(i + 1), "binary")
      
      # Set sense for the Linear Problem
      lp.control(lpmodel, sense = "max")
      
      solve(lpmodel)
    }
  }
  
  oe <- list(eff = e, 
             lambda = lambda, 
             objval = objval, 
             RTS = RTS, 
             ORIENTATION = ORIENTATION)
  
  return(oe)
  
}