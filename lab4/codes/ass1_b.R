gammh<-function (Nsim, a, b)
## Nism is the number if iterations and a and b are the parameters of the Gamma distribution 

{
  mu <- a/b
  ## the mean of the proposal distribution
  sig <- sqrt(a/(b * b))
  ## the standard deviation of the proposal distribution
  vec <- vector("numeric", Nsim)
  x <- a/b
  vec[1] <- x
  ## initializing
  accepts = 0
  for (i in 2:Nsim) {
    can <- rnorm(1, mu, sig)
    aprob <- min(1, (dgamma(can, a, b)/dgamma(x,
                                              a, b))/(dnorm(can, mu, sig)/dnorm(x,mu, sig)))
    ## computing the Metropolis-Hastings ratio
    ## here we use Metropolis-Hastings independence sampling algorithm
    ## so g(x|y) = g(x)

    u <- runif(1)
    ## ramdom probability that accepts the Xstar
    if (u < aprob){
      x <- can
    }
    vec[i] <- x
    accepts <- accepts+1
  }
  accepts
  return(vec)
}

output <- function(iteration=2000,a,b){
  vec <- gammh(iteration,a,b)
  h<-hist(vec,freq=FALSE,las=1,col=8,xlab="x", 
          main="Histogram(a = 0.1, b = 2) with Gamma Curve") 
  curve(dgamma(x, shape=0.1, scale=2), add=TRUE, col=2, lwd=2) 
  return (vec)
}
h = output(2000,0.1,2)
plot(h,las=1,col=4,xlab='iteration',ylab='chain',, 
     main="Chain(a = 0.1, b = 2) vs Iteration")
