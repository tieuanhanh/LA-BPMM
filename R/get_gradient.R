get_gradient <- function(phi, inv.cov.mat){
  # return the gradient vector of log likelihood of model given eta and data
  # phi: vector of fixed effects and random effects
  # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
  - colSums(t*y - t*c(exp(t%*%phi))) + inv.cov.mat%*%phi
}
