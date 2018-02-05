###Generate a mixture

require(grDevices) # for colours
require(mvtnorm)
x.points <-seq(-3,3,length.out=100)
y.points <-x.points
z<-matrix(0,nrow=100,ncol=100)
mu1<-c(1,1)                         ###centers of bivariate distributions
mu2<-c(-1,-1)
sigma <-matrix(c(0.5,0,0,0.1),nrow=2)   ###covariance matrix(fixed)
w1<-0.7                             ###mixture weights
w2<-1-w1
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <-rnorm(1,0,0.01)+w1*dmvnorm(c(x.points[i],y.points[j]),mean=mu1,sigma=sigma)+w2*dmvnorm(c(x.points[i],y.points[j]),mean=mu2,sigma=sigma)
  }
}
filled.contour(x.points,y.points,z,color = heat.colors)
#plot(x.points,z)

###Other simulation
x <- rnorm(200,-1,0.5)
y <- rnorm(200,1,0.5)
s <- subplot(
  plot_ly(x = x, type = "histogram"),
  plotly_empty(),
  plot_ly(x = x, y = y, type = "histogram2dcontour"),
  plot_ly(y = y, type = "histogram"),
  nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
  shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
)
p <- layout(s, showlegend = FALSE)
p


library(coda)
library(mvtnorm)

# Laplace approximation
laplace_approx <- function(model, inits, no_samples, ...) {
  fit <- optim(inits, model, control = list(fnscale = -1), hessian = TRUE, ...)
  param_mean <- fit$par
  print(param_mean)
  param_cov_mat <- solve(-fit$hessian)
  mcmc(rmvnorm(no_samples, param_mean, param_cov_mat))
  #param_cov_mat
}

sigma<-1
ndata<-100000
u<-runif(ndata)<w1
y<-u*rnorm(ndata, mean = -1, sd = sigma)+(1-u)*rnorm(ndata, mean = 1, sd = sigma)



model <- function(p, y) {
  #log_lik <- sum(dnorm(y, p["mu"], p["sigma"], log = T))
  #log_post <- log_lik + dnorm(p["mu"], 0, 100, log = T) + dlnorm(p["sigma"],0, 4, log = T)
  #log_post
  log_lik<-sum(w1*dnorm(y, p["mu1"], sigma)+(1-w1)*dnorm(y, p["mu2"], sigma, log = T))
  log_post <- log_lik + dnorm(p["mu1"], 0, 100, log = T) + dnorm(p["mu2"], 0, 100, log = T) 
  log_post
}

inits <- c(mu1 = 0, mu2 = 1.13)
samples <- laplace_approx(model, inits, 1000, y = y)
densityplot(samples)


