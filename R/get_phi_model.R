get_phi_model <- function(norm.type = "2", grad.tol = 1e-5, 
                          phi.start = rep(1,numFixEff+numUnit), 
                          max.iter = 1e5, inv.cov.mat){
  
  # return the MLE of phi given data and eta
  # norm.type is the type of distance in function "norm"
  # phi.start is our starting guess of phi. Default: vector of ones
  # grad_tol is convergence tolerance of gradient
  # max.iter is maximum number of iterations
  
  invcovmat <- inv.cov.mat
  
  # set phi be our current guess
  x <- phi.start
  
  for (i in 1:max.iter){
    # check local gradient
    j <- get_gradient(phi = x, inv.cov.mat = invcovmat)
    grad.norm <- norm(j, type = norm.type)
    if (grad.norm < grad.tol){
      return(x)
    }
    H <- get_hessian(phi = x, inv.cov.mat = invcovmat)
    inv.H <- chol2inv(chol(H))
    x <- x - inv.H%*%j
  }
  
  return(x)
}
