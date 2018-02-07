###Generate a mixture

#plot(x.points,z)

###Other simulation
# x <- rnorm(200,-1,0.5)
# y <- rnorm(200,1,0.5)
# s <- subplot(
#   plot_ly(x = x, type = "histogram"),
#   plotly_empty(),
#   plot_ly(x = x, y = y, type = "histogram2dcontour"),
#   plot_ly(y = y, type = "histogram"),
#   nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
#   shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
# )
# p <- layout(s, showlegend = FALSE)
# p


library(coda)
library(mvtnorm)
library(lattice)

model <- function(p, y) {
  #log_lik<-sum(w1*dnorm(y, p["mu1"], sigma)+(1-w1)*dnorm(y, p["mu2"], sigma, log = T))
  log_lik<-sum(log(w1*dnorm(y, p["mu1"], sigma)+(1-w1)*dnorm(y, p["mu2"], sigma)))
  log_post <- log_lik + dnorm(p["mu1"], 1, 100, log = T) + dnorm(p["mu2"], 1, 100, log = T) 
  log_post
}

# Laplace approximation
laplace_approx <- function(model, inits = c(mu1 = 0, mu2 = 1.13), no_samples, ...) {
  fit <- optim(inits, model, control = list(fnscale = -1), hessian = TRUE, ...)
  param_mean <- fit$par
  #print(param_mean)
  param_cov_mat <- solve(-fit$hessian)
  #print(round(param_cov_mat,8))
  return(list("PosteriorMean"=param_mean,"CovMat"=param_cov_mat,"samp"=mcmc(rmvnorm(no_samples, param_mean, param_cov_mat))))
}

#samples <- laplace_approx(model, inits, 100000, y = y)
#densityplot(samples$samp)
