Epanechnikov=function(x){
  
  Ke <- function(tt)	{
    sqrt(2)*(1- (tt^2)/5)
  }
  
  mytt <- seq(from=-sqrt(5), to=sqrt(5), length=length(x))
  
  sf <- sum(Ke(mytt))
  
  return(Ke(mytt)/ sf)
  
}