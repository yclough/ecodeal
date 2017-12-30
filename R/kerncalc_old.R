#' Exponential kernel computation, old version
#' 
#' An older version of kerncalc, should not be used. Kept for documentation.
#' @param beta the mean of the exponential distribution
#' @param cell.size the size of a lattice's cell in meters
#' @param maxSize The maximum alowable size of the resulting matrix.
#' @param decaycut the quantile used to cut the tails of the kernel



kerncalc_old <- function(beta,
                         cell.size,
                         maxSize = Inf,
                         decaycut=0.99) 
{
  radius = min(maxSize,round(qexp(decaycut,cell.size / beta)))
  
  dmat <- matrix(0,nrow = (radius * 2) + 1,ncol = (radius * 2) + 1)
  reachable <- matrix(0,nrow = (radius * 2) + 1,ncol = (radius * 2) + 1)
  decay <- matrix(0,nrow = (radius * 2) + 1,ncol = (radius * 2) + 1)
  
  for (ii in 1:nrow(dmat))
  {
    for (jj in 1:ncol(dmat))
    {
      reachable[ii,jj] <-
        ifelse(sqrt((ii - (radius + 1)) ^ 2 + (jj - (radius + 1)) ^ 2) > radius,0,1)
      decay[ii,jj] <-
        exp(- (cell.size/beta) * (sqrt((ii - (radius + 1)) ^ 2 + 
                                         (jj - (radius + 1)) ^ 2)))
    }
  }
  decay <- decay * reachable
  decay = decay / sum(decay)
  
  return(list(decay=decay,decaycut=decaycut))
}
