#' Spatio-Temporal DEA efficiency
#' 
#' Estimates a ST-DEA frontier and calculates efficiency measures a la Farrell.
#' @param X Inputs of firms to be evaluated, a K x m matrix of observations of K firms with m inputs (firm x input). In case TRANSPOSE=TRUE the input matrix is transposed to input x firm.
#' @param Y Outputs of firms to be evaluated, a K x n matrix of observations of K firms with n outputs (firm x input). In case TRANSPOSE=TRUE the output matrix is transposed to output x firm.
#' @param RTS Text string or a number defining the underlying DEA technology / returns to scale assumption.
#' @param ORIENTATION Input efficiency "in" (1), output efficiency "out" (2), and graph efficiency "graph" (3). For use with DIRECT , an additional option is "in-out" (0).
#' @param stp the STEP..
#' @return The results are returned in a Farrell object with the following components. The last three components in the list are only part of the object when SLACK=TRUE.
#' 
#' @export
stDEA <- function(X, Y, RTS = "vrs", ORIENTATION = "out", stp = 0.01){
  
  # imported.data <- read.csv(file = file.choose(), header = TRUE, sep = ";")
  # x <- with(imported.data, cbind(X1))
  # y <- with(imported.data, cbind(O1, O2))
  
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
  
  numberOfDMUs <- dim(X)[1]
  numberOfOutputs <- dim(Y)[2]
  numberOfInputs <- dim(X)[2]
  
  if(minimumDMUs(numberOfInputs, numberOfOutputs) > numberOfDMUs){
    stop("Not enough Decision Making Units, you should have at least : ", minimumDMUs)
  }
  
  e <- dea(X, Y, RTS, ORIENTATION)
  efficiency.dea <- e$eff
  prs <- peers(e)
  lmbd <- as.data.frame(lambda(e))
  
  alphaMatrix <- alpha.Matrix(e)
  deltaMatrix <- delta.Matrix(alphaMatrix)
  alphaMax <- alpha.Max(alphaMatrix)
  deltaMin <- delta.min(deltaMatrix)
  
  # Define the constraintTypes
  constraintTypes <- matrix(data = ">=", nrow = (numberOfOutputs + 1), ncol = 1)
  
  # Define the rhs
  rightHandSide <- matrix(data = 0, nrow = numberOfOutputs, ncol = 1)
  
  st.peer <- NULL
  #technical.efficiency <- matrix(data = NA, nrow = (1 / 0.01 + 1), ncol = numberOfOutputs)
  eff.dea <- eff(e)[numberOfDMUs]
  eff.stdea <- NULL
  technical.eff.dea <- 1 / eff.dea
  technical.eff.stdea <- NULL
  slacks <- matrix(data = NA, nrow = (1 / 0.01 + 1), ncol = numberOfOutputs)
  proj.val <- matrix(data = NA, nrow = (1 / 0.01 + 1), ncol = numberOfOutputs)
  
  # ST-DEA
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
    
    st.peer <- c(st.peer, which.max(get.variables(lpmodel)[2:i]))
    
    row.Number <- Wsp * 100 + 1
    for(ii in 1:numberOfOutputs){
      technical.efficiency[row.Number, ii] <- Y[(which.max(get.variables(lpmodel)[2:i])), ii] / Y[i, ii]
    }
    
    eff.stdea <- c(eff.stdea, min(technical.efficiency[row.Number,]))
    
    technical.eff.stdea <- c(technical.eff.stdea, 1 / eff.stdea[row.Number])
    
    for(ii in 1:numberOfOutputs){
      slacks[row.Number, ii] <- Y[(which.max(get.variables(lpmodel)[2:numberOfDMUs])), ii] - eff.stdea[row.Number]  * Y[numberOfDMUs, ii]
      proj.val[row.Number, ii] <- eff.stdea[numberOfDMUs] * Y[numberOfDMUs, ii]
    }
    
    rm(lpmodel)
  }
  View(eff.dea)
  View(eff.stdea)
  View(technical.eff.dea)
  View(technical.eff.stdea)
  View(slacks)
  View(proj.val)
  
  oe <- list(effDEA = eff.dea,
             effstDEA = eff.stdea,
             lambda = lambda,
             objval = objval,
             RTS = RTS,
             ORIENTATION = ORIENTATION)
  
  return(oe)
  
}