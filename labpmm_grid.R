labpmm_grid <- function(y, t, numUnit, obsPerUnit, numFixEff,
                        gammaA = 0.01, gammaB = 0.01, fixEffVar = 10^4,
                        grid_n = 10){
  library(Matrix)
  # y: the response variable (count data)
  # t: the predictive variables
  # numUnit: number of groups/units, aka number of random intercepts 
  # in case of one random effect intercept for each group/unit
  # obsPerUnit: number of observations per group/unit
  # numFixEff: number of fixed effect variables
  # gammaA, gammaB: parameters of the prior precision eta ~ Gamma(a,b)
  # fixEffVar: variance of the prior of fixed effect B ~ N(0, fixEffVar*I)
  
  # Design matrix for the random intercept
  # z1: column vector of ones
  z1 = matrix(rep(1,obsPerUnit), ncol=1)
  # z: design matrix of random effect intercept b0
  z = kronecker(diag(numUnit), z1)
  ones = rep(1, length(y))
  t = cbind(ones, t,z)
  t = as.matrix(t)
  
  ## Calculate likelihood of data given eta
  
  get_minus_log_llh <- function(phi, inv.cov.mat){
    # return minus log likelihood of model given eta and data
    # phi: vector of fixed effects and random effects
    # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
    - sum(t %*% phi * y - exp(t%*%phi)) + 1/2 * phi %*% inv.cov.mat %*% phi
  }
  
  ## Gradient and hessian of f
  
  get_gradient <- function(phi, inv.cov.mat){
    # return the gradient vector of log likelihood of model given eta and data
    # phi: vector of fixed effects and random effects
    # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
    - colSums(t*y - t*c(exp(t%*%phi))) + inv.cov.mat%*%phi
  }
  
  get_hessian <- function(phi, inv.cov.mat){
    # return the hessian matrix of log likelihood of model given eta and data
    # phi: vector of fixed effects and random effects
    # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
    t(t)%*%(t*c(exp(t%*%phi))) + inv.cov.mat
  }
  
  ##
  
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
  
  get_func_for_v <- function(v, gamma.a, gamma.b, fix.eff.var){
    ## Return function of posterior of eta 
    .cov.mat = diag(rep(c(fix.eff.var, 1/exp(v)), times = c(numFixEff,numUnit)))
    .inv.cov.mat <- chol2inv(chol(.cov.mat))
    .phi.mode <- get_phi_model(inv.cov.mat = .inv.cov.mat)
    .H <- get_hessian(phi = .phi.mode, inv.cov.mat = .inv.cov.mat)
    .la.var <- chol2inv(chol(.H))
    funcval <- get_minus_log_llh(c(.phi.mode), inv.cov.mat = .inv.cov.mat) + 
      gamma.b * exp(v) - (numUnit/2+gamma.a)*v -1/2* log(det(.la.var))
    return(funcval)
  }
  
  vMode <- optimize(get_func_for_v, c(-20, 20), tol = 0.0001, 
                    gamma.a = gammaA, gamma.b = gammaB, 
                    fix.eff.var = fixEffVar)
  vMode <- vMode$minimum
  
  ### Approach 2
  ### Laplace Approximation at vMode
  # its corresponding covariance matrix
  # covMatMode = diag(rep(c(fixEffVar, 1/exp(vMode)), times = c(numFixEff,numUnit)))
  # invCovMatMode <- chol2inv(chol(covMatMode))
  # # estimation of fixed effects and random effects
  # phiModel <- get_phi_model(phi.start = rep(1,numFixEff+numUnit),
  #                           inv.cov.mat = invCovMatMode)
  # fixEffectModel <- phiModel[1:numFixEff]
  # randomEffectModel <- phiModel[numFixEff+1:numUnit]
  # 
  # hessianMode <- get_hessian(phiModel, inv.cov.mat = invCovMatMode)
  # la.var.model <- chol2inv(chol(hessianMode))
  # 
  # diag = la.var.model[row(la.var.model) == col(la.var.model)]
  # sd_fixEffect = sqrt(diag[1:numFixEff])
  
  #### Approach 1
  i = 0
  v_grid = seq(vMode - grid_n/2*0.15, vMode + grid_n/2*0.15, 0.15)
  log_post_v = rep(0, length(v_grid))
  phi.mode = matrix(nrow =  length(v_grid), ncol = numFixEff)
  phi.var = matrix(nrow =  length(v_grid), ncol = numFixEff)
  
  for (v in v_grid){
    i = i + 1
    log_post_v[i] = -get_func_for_v(v, gamma.a = gammaA, 
                                    gamma.b = gammaB, 
                                    fix.eff.var = fixEffVar)
    .cov.mat = diag(rep(c(fixEffVar, 1/exp(v)), times = c(numFixEff,numUnit)))
    .inv.cov.mat <- chol2inv(chol(.cov.mat))
    phi.vector <- get_phi_model(inv.cov.mat = .inv.cov.mat)
    .H <- get_hessian(phi = phi.vector, inv.cov.mat = .inv.cov.mat)
    la.var <- chol2inv(chol(.H))
    for (j in 1:numFixEff){
      phi.mode[i,j] = phi.vector[j]
      phi.var[i,j] = la.var[j,j]
    }
  }
  
  log_post_v_scale = log_post_v - max(log_post_v)
  post_v_scale = exp(log_post_v_scale)
  post_v_norm = post_v_scale/sum(post_v_scale)
  
  int.est <- colSums(phi.mode * post_v_norm)
  est.var <- colSums(phi.var * post_v_norm)
  
  # Posterior mean of the variance
  var.est <- sum(1/exp(v_grid) * post_v_norm)
  m2_var = sum((1/exp(v_grid)^2)*post_v_norm)
  sd_var = sqrt(m2_var - var.est^2)
  
  ### get the CI of variance
  mean_v = sum(v_grid*post_v_norm)
  m2_v = sum(v_grid^2*post_v_norm)
  var_v = m2_v - mean_v^2
  sd_v = sqrt(var_v)
  lb_95_v = mean_v - 1.96*sd_v
  ub_95_v = mean_v + 1.96*sd_v
  ub_95_var = 1/exp(lb_95_v)
  lb_95_var = 1/exp(ub_95_v)
  ci_95_var = c(lb_95_var, ub_95_var)
  
  lb_90_v = mean_v - 1.645*sd_v
  ub_90_v = mean_v + 1.645*sd_v
  ub_90_var = 1/exp(lb_90_v)
  lb_90_var = 1/exp(ub_90_v)
  ci_90_var = c(lb_90_var, ub_90_var)
  
  results <- list(fix_effect = int.est,
                  sd_fixEffect = sqrt(est.var), 
                  variance = var.est,
                  sd_var = sd_var,
                  ci_95_var = ci_95_var, 
                  ci_90_var = ci_90_var,
                  vMode = vMode, sd_v = sd_v, mean_v = mean_v)
  return(results)
}
