#' LMED
#'
#' @param x Covariate
#' @param y response
#' @param loc_age location of the age in the covariate matrix
#' @param lambda regularizer for the stability of the algorithm
#'
#' @return geometric median of planar shape
#' @export

LMED <- function(x,y,age_loc,lambda){

  age <- x[,age_loc]
  N <- dim(y)[1] # total number of observations
  S <- dim(y)[2] # total number of locations
  P <- dim(x)[2]+1 # dimension of the design matrix (intercept included)


  ite <- 100 # maximum number of iterations
  Tol <- 1e-6


  delta <- Delta_modi(x,y,age_loc,lambda)
  brain_age <- age+delta

  x_delta <- data.frame(x)
  x_delta[,age_loc] <- brain_age # Covariate matrix plugging age with brain age

  # Initializing Step

  psi <- rnorm(N,1,0.01) # Initialize the random coefficient
  Psi <- diag(psi)

  temp <- lm(y~x)
  coef <- coef(temp) # Beta_0
  pred <- predict(temp)
  pred.delta <- predict(temp,newdata=x_delta)
  G <- pred.delta-pred # initial value of G

  # Initialize covariance components
  Omega <- t(Psi%*%G)%*%(Psi%*%G)/N
  E <- resid(temp)-G
  Sigma <- t(E)%*%E/N
  V <- Omega+Sigma + diag(lambda,S)

  # Some setup for the main iterations

  one <- rep(1,S)
  I <- diag(N)

  X.kro <- kronecker(one,cbind(1,x))
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
    predict.delta <- as.matrix(cbind(1,x_delta))%*%B_hat
    predict.y <- cbind(1,x)%*%B_hat
    G.new <- predict.delta - predict.y

    # Update Psi
    psi.mean <- mean(psi)
    psi.sig <- var(psi)
    psi.new <- psi.mean + psi.sig*diag(G.new%*%solve(V)%*%t(y-cbind(1,x)%*%B_hat))

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

  estimated <- cbind(1,x)%*%B_hat.new + diag(psi.new)%*%G.new

  # Estimate the SE of Reg Coeff
  Sig.hat <- colSums((y-estimated)^2)/(N-P)
  Xtx.inv <- diag(solve(t(cbind(1,x))%*%cbind(1,x)))
  Se.beta <- sapply(1:S, function(i) sqrt(Sig.hat[i]*Xtx.inv))


  # Test Statistic
  test.statistic <- B_hat.new/Se.beta
  P.value <- pt(abs(test.statistic),df=N-P,lower.tail=FALSE)


  result <- list("fitted"=estimated,"coef"=B_hat.new,"SE"=Se.beta,"Pval"=P.value)
  return(result)
  return(B_hat.new)
}


