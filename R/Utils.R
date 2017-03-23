### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## Weighting around the precursor
#' 
#' Weighting around the precursor window -> to be improved with larger windo etc...
#' 
#' @param x vector of mz-mzPrec
#' @param weight weight
#' @keywords internal
#' 
#' @export
.MGgetMZweight<-function(x,weight){
  if(is.na(weight)) return(rep(1,length(x)))
  we= 1-abs(weight)*x^2
  we[we<=0]=0
  we
}
