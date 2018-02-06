## EP for Gaussian mixture

## prior p(theta) = N(0,I)

EPGM = function(n,x,iter,likelihood_var = 1,w){
  ## n: number of observations
  ## x: observations
  ## iter: number of iterations
  ## likelihood_var: variance of xi
  ## w: weight for the first mixture
  
  ## initialize every ti(theta) to be N(0,I)
  mu_initial = matrix(rep(0,2*n),2,n)
  Sigma_initial = list(0)
  for(j in 1:n){
    Sigma_initial[[j]] = matrix(c(1,0,0,1),2,2)
  }
  
  ## the initial q(theta)
  mu = matrix(c(0,0),2,1)
  Sigma = (1/(1+n))*matrix(c(1,0,0,1),2,2)
  
  mu_all = mu
  Sigma_all = Sigma
  
  ## Update every ti until convergence
  for(ind in 1:iter){
    
    for(i in 1:n){
      
      ## cavity distribution
      mu_tilt = as.matrix(mu_initial[,i],2,1)
      Sigma_tilt = Sigma_initial[[i]]
      
      Sigma__i = solve(solve(Sigma) - solve(Sigma_tilt))
      mu__i = Sigma__i %*% (solve(Sigma) %*% mu - solve(Sigma_tilt) %*% mu_tilt)
      
      ## recompute the approximation combining the real likelihood
      Sigma_star = solve(solve(Sigma__i) + matrix(c(likelihood_var,0,0,0),2,2))
      mu_star = Sigma_star %*% (solve(Sigma__i) %*% mu__i + matrix(c(likelihood_var,0,0,0),2,2) %*% matrix(c(x[i],0),2,1))
      
      Sigma_hat = solve(solve(Sigma__i) + matrix(c(0,0,0,likelihood_var),2,2))
      mu_hat = Sigma_star %*% (solve(Sigma__i) %*% mu__i + matrix(c(0,0,0,likelihood_var),2,2) %*% matrix(c(0,x[i]),2,1))
      
      mu_prime = w * mu_star + (1-w) * mu_hat
      Sigma_prime = w * Sigma_star + (1-w) * Sigma_hat + w*(1-w)* (mu_star - mu_hat) %*% t(mu_star - mu_hat)
      
      mu = mu_prime
      Sigma = Sigma_prime
      
      ## update ti
      Sigma_new = solve(solve(Sigma_prime) - solve(Sigma__i))
      mu_new = Sigma_new %*% (solve(Sigma_prime) %*% mu_prime - solve(Sigma__i) %*% mu__i)
      
      Sigma_initial[[i]] = Sigma_new
      mu_initial[,i] = mu_new
    }
    mu_all = cbind(mu_all,mu)
    Sigma_all = cbind(Sigma_all,Sigma)
    
  }
  return(list(mu = mu_all,Sigma = Sigma_all))
}


## sample the observations
n = 100
mu_prior = c(0,0)
var_prior = c(1,1)
theta1 = rnorm(1,mu_prior[1],var_prior[1])
theta2 = rnorm(1,mu_prior[2],var_prior[2])
w = sample(c(1,0),n,replace = TRUE)
theta = c(theta1,theta2)
x = w * rnorm(1,theta1,1) + (1-w) * rnorm(1,theta2,1)

## initialise parameters
n = 100
x = x
iter = 100
likelihood_var = 1
w = 0.5

EP1 = EPGM(n = 100,x = x,iter = 400,likelihood_var = 1,w = 0.5)
EP1

