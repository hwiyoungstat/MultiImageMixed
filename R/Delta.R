#' Delta_modi
#'
#' @param x Covariate
#' @param y response
#' @param loc_age location of the age in the covariate matrix
#' @param lambda regularizer for the stability of the algorithm
#'
#' @return geometric median of planar shape
#' @export

Delta_modi <- function(x,y,loc_age,lambda){

  z <- as.matrix(cbind(1,y,x[,-loc_age]))
  age <- as.matrix(x[,loc_age]) # used for regression 2 (age as response)
  age1 <- as.matrix(cbind(1,x[,loc_age])) # used for regression 1 (age as covriate)

  n <- dim(y)[1]
  s <- dim(z)[2]

  temp1 <- z%*%solve(t(z)%*%z + diag(lambda,s))%*%t(z)%*%age-age
  temp2 <- temp1-age1%*%solve(t(age1)%*%age1)%*%t(age1)%*%temp1

  return(temp2)
}




