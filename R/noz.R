#' No zero helper function
#' 
#' Replace negative values by zeros and rescale such that the same sum is preserved
#' @param x a matrix or a vector
#' @return same as x, that is a matrix or a vector
noz <- function(x){
  s<-sum(x)
  x[x<0]<-0
  x<-x/sum(x)*s
}