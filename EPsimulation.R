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

# Laplace approximation
laplace_approx <- function(model, inits, no_samples, ...) {
  fit <- optim(inits, model, control = list(fnscale = -1), hessian = TRUE, ...)
  param_mean <- fit$par
  print(param_mean)
  param_cov_mat <- solve(-fit$hessian)
  #print(round(param_cov_mat,8))
  return(list("PosteriorMean"=param_mean,"CovMat"=param_cov_mat,"samp"=mcmc(rmvnorm(no_samples, param_mean, param_cov_mat))))
}

sigma<-0.5
ndata<-100000
set.seed(332)
u<-runif(ndata)<w1
y<-u*rnorm(ndata, mean = -3, sd = sigma)+(1-u)*rnorm(ndata, mean = 3, sd = sigma)



model <- function(p, y) {
  #log_lik<-sum(w1*dnorm(y, p["mu1"], sigma)+(1-w1)*dnorm(y, p["mu2"], sigma, log = T))
  log_lik<-sum(log(w1*dnorm(y, p["mu1"], sigma)+(1-w1)*dnorm(y, p["mu2"], sigma)))
  log_post <- log_lik + dnorm(p["mu1"], 1, 100, log = T) + dnorm(p["mu2"], 1, 100, log = T) 
  log_post
}

inits <- c(mu1 = 0, mu2 = 1.13)
samples <- laplace_approx(model, inits, 100000, y = y)
densityplot(samples$samp)



require(grDevices) # for colours
require(mvtnorm)
x.points <-seq(-5,5,length.out=100)
y.points <-x.points
z.post<-matrix(0,nrow=100,ncol=100)
z.lap<-z.post
###Posterior distribution
mu1<-c(3,-3)                         ###centers of bivariate distributions
mu2<-c(-3,3)
sigma.mat <-matrix(c(sigma,0,0,sigma),nrow=2)   ###covariance matrix(fixed)
w1<-0.5                             ###mixture weights
w2<-1-w1
for (i in 1:100) {
  for (j in 1:100) {
    z.post[i,j] <-w1*dmvnorm(c(x.points[i],y.points[j]),mean=mu1,sigma=sigma.mat)+w2*dmvnorm(c(x.points[i],y.points[j]),mean=mu2,sigma=sigma.mat)
  }
}
filled.contour(x.points,y.points,z.post,color = terrain.colors,asp=1)

for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <-dmvnorm(c(x.points[i],y.points[j]),mean=c(-3,3),sigma=sigma.mat)
    z.lap[i,j] <-dmvnorm(c(x.points[i],y.points[j]),mean=samples$PosteriorMean,sigma=sigma.mat)
  }
}
contour(x.points,y.points,z,add=TRUE,col=2)
contour(x.points,y.points,z.lap,add=TRUE)
