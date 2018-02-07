###Contour plots
require(grDevices) # for colours
require(mvtnorm)

###Generating the data
n.iter<-100
sigma<-0.5
ndata<-100
w1<-0.5                             ###mixture weights
w2<-1-w1
u<-runif(ndata)<w1
y<-u*rnorm(ndata, mean = -3, sd = sigma)+(1-u)*rnorm(ndata, mean = 3, sd = sigma)
hist(y)

###Generating the plot grid
x.points <-seq(-5,5,length.out=100)
y.points <-x.points
z.post<-matrix(0,nrow=100,ncol=100)
z.lap<-z.post
z.vb<-z.post
z.ep<-z.post


###Generating the approximations to the posterior
estLAP<-laplace_approx(model, inits = c(mu1 = 0, mu2 = 1.13), 100, y = y)
estVB<-VBMixture(data = y, sigma = sigma,  eps = 1e-3)
estEP<-EPGM(n = length(y),x = y,iter = n.iter,likelihood_var = sigma,w = 0.5)

estVB.mean<-estVB[1,]
estVB.CovMat<-matrix(c(estVB[2,1],0,0,estVB[2,2]),2,2)
estEP.mean<-estEP$mu[,dim(estEP$mu)[2]]
#estEP.CovMat<-round(matrix(c(estEP$Sigma[,(dim(estEP$Sigma)[2]-1):dim(estEP$Sigma)[2]]),2,2))
estEP.CovMat<-matrix(c(estEP$Sigma[,(dim(estEP$Sigma)[2]-1):dim(estEP$Sigma)[2]]),2,2)
#estEP.CovMat<-matrix(c(estEP$Sigma[,3:4]),2,2)


model <- function(p, y) {
  #log_lik<-sum(w1*dnorm(y, p["mu1"], sigma)+(1-w1)*dnorm(y, p["mu2"], sigma, log = T))
  log_lik<-sum(log(w1*dnorm(y, p["mu1"], sigma)+(1-w1)*dnorm(y, p["mu2"], sigma)))
  log_post <- log_lik + dnorm(p["mu1"], 1, 100, log = T) + dnorm(p["mu2"], 1, 100, log = T) 
  return("logpost"=log_post)
}

for (i in 1:100) {
  for (j in 1:100) {
    z.post[i,j]<-model(c("mu1"=x.points[i],"mu2"=y.points[j]),y)
    z.lap[i,j] <-dmvnorm(c(x.points[i],y.points[j]),mean=estLAP$PosteriorMean,sigma=estLAP$CovMat)
    z.vb[i,j] <-dmvnorm(c(x.points[i],y.points[j]),mean=estVB.mean,sigma=estVB.CovMat)
    z.ep[i,j] <-dmvnorm(c(x.points[i],y.points[j]),mean=estEP.mean,sigma=estEP.CovMat) 
  }
}
contour(x.points,y.points,exp(z.post)/sum(exp(z.post)))#,color=terrain.colors,colorkey=FALSE)
contour(x.points,y.points,z.lap,add=TRUE, col=2)
contour(x.points,y.points,z.vb,add=TRUE, col=5)
contour(x.points,y.points,z.ep,add=TRUE, col=3)

