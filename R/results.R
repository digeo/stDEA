oe <- list(eff = e, lambda = lambda, objval = objval, RTS = RTS, 
           primal = primal, dual = dual, ux = ux, vy = vy, gamma = gamma, 
           ORIENTATION = ORIENTATION, TRANSPOSE = TRANSPOSE, param = param)
if (!is.null(DIRECT)) {
  oe$direct <- DIRECT
}
class(oe) <- "Farrell"
if (SLACK) {
  if (TRANSPOSE) {
    X <- t(X)
    Y <- t(Y)
    if (.xyref.missing) {
      XREF <- NULL
      YREF <- NULL
    }
    else {
      XREF <- t(XREF)
      YREF <- t(YREF)
    }
  }
  sl <- slack(X, Y, oe, XREF, YREF, FRONT.IDX, LP = LP)
  oe$slack <- sl$slack
  oe$sum <- sl$sum
  oe$sx <- sl$sx
  oe$sy <- sl$sy
  oe$lambda <- sl$lambda
  if (LP) {
    print("slack fra slack:")
    print(sl$slack)
    print("slack efter slack:")
    print(oe$slack)
  }
}
return(oe)