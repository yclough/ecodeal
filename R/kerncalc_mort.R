#' Exponential kernel computation, with mortality
#' 
#' Exponential kernel computation to be used as an input of the latfordisp function.
#' @param beta the mean of the exponential distribution
#' @param cell.size the size of a lattice's cell in meters
#' @param maxSize The maximum alowable size of the resulting matrix.
#' @param decaycut the quantile used to cut the tails of the kernel
#' @param surv.per.m proportion of individuals surviving per meter travelled

kerncalc_mort <- function(beta, cell.size, maxSize = Inf, decaycut=0.99,surv.per.m=0.9)
{
  radius = min(maxSize, round(qexp(decaycut,cell.size / beta)))
  mat_nrow = (radius * 2) + 1
  mat_ncol = (radius * 2) + 1
  ij <- cbind(1:mat_nrow, rep(1:mat_ncol, each= mat_nrow))
  reachable <- sqrt((ij[,1] - (radius + 1)) ^ 2 + (ij[,2] - (radius + 1)) ^ 2) <= radius
  decay <- exp(- (cell.size/beta) * (sqrt((ij[,1]  - (radius + 1)) ^ 2 + (ij[,2]  - (radius + 1)) ^ 2)))
  decay <- decay * reachable
  decay = decay / sum(decay)
  decay <- matrix(decay, mat_nrow, mat_ncol)
  mortality <- surv.per.m^( cell.size *(sqrt((ij[,1]  - (radius + 1)) ^ 2 + (ij[,2]  - (radius + 1)) ^ 2)))
  return(list(decay=decay*mortality,decaycut=decaycut))
}
