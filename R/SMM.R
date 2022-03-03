#' SMED
#'
#' @param x Covariate
#' @param y response
#' @param loc_age location of the age in the covariate matrix
#' @param lambda regularizer for the stability of the algorithm
#' @return geometric median of planar shape
#' @export

SMED <- function(x,y,age_loc,lambda){

  age <- x[,age_loc]
  N <- dim(y)[1] # total number of observations
  S <- dim(y)[2] # total number of locations


  ite <- 100 # maximum number of iterations
  Tol <- 1e-6


  delta <- Delta_modi(x,y,age_loc,lambda)
  brain_age <- age+delta

  min.knot <- min(brain_age,age)
  max.knot <- max(brain_age,age)

  # Setup for Spline

  BS <- ns(age,df=2,Boundary.knots=c(min.knot,max.knot))
  Delta.BS <-  ns(brain_age,knots=attr(BS,"knots"),Boundary.knots=attr(BS,"Boundary.knots"))

  Q <- dim(BS)[2] # number of knots regarding age effect (# of pars for spline)
  M <- dim(x)[2]-1 # dimension of covariate except age
  P <- Q+M+1 # number of parameters for final spline regression (including intercept)


  x.BS <- cbind(x[,-age_loc],BS)
  x.delta.BS <- data.frame(cbind(x[,-age_loc],Delta.BS))


  # Initializing Step

  psi <- rnorm(N,1,0.01) # Initialize the random coefficient
  Psi <- diag(psi)

  temp <- lm(y~x.BS)
  coef <- coef(temp) # Beta_0
  pred <- predict(temp)
  pred.delta <- predict(temp,newdata=x.delta.BS)
  G <- pred.delta-pred # initial value of G



  # Initialize covariance components
  Omega <- t(Psi%*%G)%*%(Psi%*%G)/N
  E <- resid(temp)-Psi%*%G
  Sigma <- t(E)%*%E/N
  V <- Omega+Sigma + diag(lambda,S)

  # Some setup for the main iterations
  V.inv <- solve(V)

  one <- rep(1,S)
  I <- diag(N)

  X.kro <- kronecker(one,cbind(1, x.BS))

  Y.kro <- kronecker(one,y)
  V.inv <- solve(V)
  Temp.array <- array(t(X.kro),c(P,S,N))
  Temp.mat <- sapply(1:N, function(i) Temp.array[,,i]%*%V.inv,simplify="array")
  Temp.mat <- matrix(Temp.mat,nrow=P)
  B_hat <- solve(Temp.mat%*%X.kro)%*%(Temp.mat%*%Y.kro)






  # Iteration

  for(i in 1:ite){
    print(i)
    # Updata G
    predict.delta <- as.matrix(cbind(1,x.delta.BS))%*%B_hat
    predict.y <- cbind(1,x.BS)%*%B_hat
    G.new <- predict.delta - predict.y

    # Update Psi
    psi.mean <- mean(psi)
    psi.sig <- var(psi)
    psi.new <- psi.mean + psi.sig*diag(G.new%*%solve(V)%*%t(y-cbind(1,x.BS)%*%B_hat))

    # Updata Covariance components
    # 1. Update Omege
    Psi.new <- diag(psi.new)
    Omega.new <- t(Psi.new%*%G.new)%*%(Psi.new%*%G.new)/N

    # 2. Update Sigma
    E.new <- predict.y-Psi.new%*%G.new
    Sigma.new <- t(E.new)%*%E.new/N

    # 3. Update V
    V.new <- Omega.new+Sigma.new + diag(lambda,S)

    # 4. Update Fixed Effects reg coefficients
    V.new.inv <- solve(V.new)

    Temp.mat <- sapply(1:N, function(i) Temp.array[,,i]%*%V.new.inv,simplify="array")
    Temp.mat <- matrix(Temp.mat,nrow=P)
    B_hat.new <- solve(Temp.mat%*%X.kro)%*%(Temp.mat%*%Y.kro)


    if(norm(B_hat.new-B_hat,type="F") <Tol) break

    # Update components for the next iterations
    B_hat <- B_hat.new
    V <- V.new
    psi <- psi.new
  }

  estimated <- cbind(1,x.BS)%*%B_hat.new +  diag(psi.new)%*%G.new


  # Estimate the SE of Reg Coeff
  Sig.hat <- colSums((y-cbind(1,x.BS)%*%B_hat.new)^2)/(N-P)
  Xtx.inv <- diag(solve(t(cbind(1,x.BS))%*%cbind(1,x.BS)))
  Se.beta <- sapply(1:S, function(i) sqrt(Sig.hat[i]*Xtx.inv))


  # Test Statistic
  test.statistic <- B_hat.new/Se.beta
  P.value <- pt(abs(test.statistic),df=N-P,lower.tail=FALSE)


  result <- list("fitted"=estimated,"coef"=B_hat.new,"SE"=Se.beta,"Pval"=P.value)


  return(result)
}
