# theta <- c(-1,1)
# sigma <- 0.25
# n <- 100
# ind <- sample(c(1,2), n, replace = TRUE)
# data <- rnorm(n, sd = sigma) + theta[ind] 



VBMixture<-function(data, sigma,  eps = 1e-3){
  
  #prior 2-dim standard normal
  Q = matrix(0:1, nrow = 2, ncol = 2) 
  colnames(Q)<-c("1st component","2nd component")
  rownames(Q)<-c("PosteriorMean","PosteriorVar")
  z = sample(1:2, length(data), rep=T) #initialize labels
  Q_old = matrix(0, nrow = 2, ncol = 2)

  while(max(abs(Q - Q_old))>eps){
    
    Q_old = Q
    
      Q[2,1] =  sigma^2/(sum(z==1) + sigma^2) #var of Q1
      Q[1,1] =  sum(z==1)/(sigma^2+sum(z==1))*mean(data[z==1]) #mean of Q1

      Q[2,2] =  sigma^2/(sum(z==2) + sigma^2) #var of Q2
      Q[1,2] =  sum(z==2)/(sigma^2+sum(z==2))*mean(data[z==2]) #mean of Q2
    
    #update labels
    #cluster 2 if density of cluster 1 lower than that of cluster 2
    z = as.numeric(dnorm(data, Q[1,1], sqrt(Q[2,1]), log = T) < dnorm(data, Q[1,2], sqrt(Q[2,2]), log = T)) + 1 

  }
  return(Q)
  
}


