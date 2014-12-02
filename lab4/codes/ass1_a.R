normm<-function (Nsim, a)
## Nism is the number if iterations and the a is 
## if we choose a = 5.3, it seems the acceptance rate is
## around 0.3
## the parameter of the Proposal distribution g(.|x^t)
## here we choose uniform distribution
{
  vec <- vector("numeric", Nsim)
  ## the vector vec is used to store the simulation result
  x <- 0
  vec[1] <- x
  ## innitialzing
  accepts = 0
  for (i in 2:Nsim) {
    innov <- runif(1, -a, a)
    Xstar <- x + innov
    ## sample a candidate value Xstar 
    ## from the proposal distribution
    aprob <- min(1, dnorm(Xstar)/dnorm(x))
    ## computing the Metropolis-Hastings ration
    ## here we using Metroplolis algorithm 
    ## so g(x|y) = g(y|x)
    u <- runif(1)
    ## ramdom probability that accepts the Xstar
    if (u < aprob){
      x <- Xstar
      accepts <- accepts +1
    }
    vec[i] <- x
  }
  print(accepts/Nsim)
return (vec)
}

output <- function(iteration=2000,alpha=5.3){
  vec <- normm(iteration,alpha)
  h<-hist(vec, xlim = c(-5,5),freq=FALSE,las=1,col=8,xlab="x", 
          main="Histogram(Alpha = 5.3) with Normal Curve") 
  curve(dnorm(x, mean=0, sd=1), add=TRUE, col=2, lwd=2) 
  return (vec)
}
h = output(2000,5.3)
plot(h,las=1,col=4,xlab='iteration',ylab='chain',, 
     main="Chain(Alpha = 5.3) vs Iteration")

