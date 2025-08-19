countends.stlpp <- function(X,
                            r,
                            t){
  
  if (!inherits(X, "stlpp")) stop("X should be from class stlpp")
  
  Y <- as.lpp.stlpp(X)
  L <- domain(Y)
  n <- npoints(Y)
  tvec <- as.numeric(X$data$t)
  tbordl <- X$time[1]
  tbordr <- X$time[2]
  
  comb <- data.frame(r,t)
  out <- list()
    
    for (i in 1:n) {
      
      out[[i]] <- unlist(lapply(X=1:nrow(comb), function(j){
        cof <- (tvec[i]-comb[j,2] > tbordl) + (tvec[i]+comb[j,2] < tbordr)
        ct <- countends(L,x=Y[i],r=comb[j,1])*cof
      }))
  }
  
  if(length(out)==1) out <- unlist(out)
  names(out) <- paste0("r=",comb[,1], "," , "t=",comb[,2])
  return(out)
}
