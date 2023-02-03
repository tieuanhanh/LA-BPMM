
get_inference_fix_effect <- function(phi.model, la.var.mode, i){
  
  ### Inference for fix effects
  # phi.model: vector of estimation of fixed effects
  # eta.model: estimation of the precision of the prior of random effects
  
  beta.mode = phi.model[i]
  beta.var = la.var.mode[i,i]
  x.grid = seq(beta.mode-1,beta.mode+1,by=0.001)
  beta.posterior = dnorm(x.grid, beta.mode, sd=sqrt(beta.var))
  lw.bound.hpd = beta.mode - qnorm(0.975)*sqrt(beta.var)
  up.bound.hpd = beta.mode + qnorm(0.975)*sqrt(beta.var)
  
  print("95 CI for LA of marginal posterior beta0")
  print(lw.bound.hpd)
  print(beta.mode)
  print(up.bound.hpd)
  
  plot(x.grid, beta.posterior, col = "black", type = "l", lty = 1, xlab = "beta",
       ylab="P(beta/eta.mode,D)", main = "LA for Marginal Posterior of beta")
  abline(v=c(lw.bound.hpd, beta.mode, up.bound.hpd), col = c("blue", "red", "blue"), lty= c(3,2,3))
  legend("topright", legend=c("lower bound CI", "true beta0", "upper bound CI"), 
         col = c("blue", "red", "blue"), lty= c(3,2,3), cex = 0.6)
}

# Example

for (i in 1:numFixEff){
  get_inference_fix_effect(phi.model = model2$fix_effect, 
                           la.var.mode = model2$la.var.model, i)
}

