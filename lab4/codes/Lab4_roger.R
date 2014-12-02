
#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
#  732A38 Computational Statistics - Computer lab 4: Markov Chain Monte Carlo
#  Roger Karlsson: 2014-05-02
#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------

setwd("C:/Users/Roger/Dropbox/732A38 - Computational statistics/Labs/4 - Markov Chain Monte Carlo")

# Packages
library(MASS)

#install.packages("MCMCpack")
library(MCMCpack)

install.packages('statmod')
library(statmod)

#install.packages('coda')
library(coda)

library(xtable)


# Assignment 1: Different versions of the Metropolis-Hastings algo --------


# (a) ---------------------------------------------------------------------



normm<-function (Nsim, a) 
{
  # A vector for the values
  vec <- vector("numeric", 10)
  
  # Start value
  x <- 0
  vec[1] <- x
  
  for (i in 2:Nsim) {
    # Draw a value from the proposal distribution
    innov <- runif(1, -a, a)
    
    # Add the runif value to previous value in the chain. This is a random walk. 
    Xstar <- x + innov
    
    # Calculate the acceptance ratio
    aprob <- min(1, dnorm(Xstar)/dnorm(x))
    # A random value between 0 and 1 to compare with
    u <- runif(1)
    
    # x is updated with the proposal value if the acceptance ratio is higher than the uniform value. 
    if (u < aprob) 
      x <- Xstar
    vec[i] <- x
  }
  
  # Return the simulated values
  vec
  
  
}

run1 <- normm(2000, .1)
run2 <- normm(2000, 1)
run3 <- normm(2000, 200)
run4 <- normm(2000, 5)


par(mfrow=c(2,2))
plot(run1, type= 'l')
plot(run2, type= 'l')
plot(run3, type= 'l')
plot(run4, type= 'l')


hist(run1)
hist(run2)
hist(run3)
hist(run4)

par(mfrow=c(1,1))


# Acceptance rate
length(unique(run1))/2000
length(unique(run2))/2000
length(unique(run3))/2000
length(unique(run4))/2000



# (b) ---------------------------------------------------------------------


gammh<-function (Nsim, a, b) 
{
  # Calculating input for the proposal distribution
  mu <- a/b  
  sig <- sqrt(a/(b * b))
  
  # A vector for storing
  vec <- vector("numeric", Nsim)
  
  # Start value
  x <- a/b
  vec[1] <- x
  
  for (i in 2:Nsim) {
    
    # Generate candidate from the normal distribution
    can <- rnorm(1, mu, sig)
    
    # Calculates the acceptance probability
    aprob <- min(1, (dgamma(can, a, b)/dgamma(x, 
                                              a, b))/(dnorm(can, mu, sig)/dnorm(x,mu, sig)))
    
    # A random value between 0 and 1
    u <- runif(1)
    # Accept the candidate if acceptance prob greater than u
    if (u < aprob) 
      x <- can
    vec[i] <- x
  }
  
  # Return the vector
  vec
}


run1 <- gammh(Nsim=10000, a=0.1, b=0.01)
run2 <- gammh(Nsim=10000, a=2, b=2)

par(mfrow= c(1,2))
plot(run1, type='l')
plot(run2, type='l')

hist(run1)
hist(run2)
par(mfrow= c(1,1))


length(unique(run1))/10000
length(unique(run2))/10000



# Assignment 2: The Gibbs sampling algorithm for the one-way rando --------

# DATA
data <- read.table("ass2.txt")

# Init values
a <- rep(0.001, 3)
b <- rep(0.001, 3)
mu0 <- 0

hGibbs <- function(burn, niter, data, a, b, mu0){
  
  k <- length(table(data$Group))
  Ybar <- by(data$Y, data$Group, mean)
  n <- by(data$Y, data$Group, length)
  ntot <- nrow(data)
  
  # Vectors and matrices for storing
  theta <- matrix(0, ncol = k, nrow = (burn+niter))
  mu <- vector('double', (burn+niter))
  sigma2 <- vector('double', (burn+niter))
  tau2 <- vector('double', (burn+niter))
  sigma2u <- vector('double', (burn+niter))
  
  sigma2[1] <- tau2[1] <- sigma2u[1] <- 1000
  mu[1] <- mean(Ybar)
  theta[1,] <- Ybar
  
  for(i in 2:(burn+niter)){
    
    # THETA
    theta[i,] <- mvrnorm(n = 1, mu= as.vector(
      (sigma2[(i-1)] / (sigma2[(i-1)] + n*tau2[(i-1)]))*mu[(i-1)] + 
        ((n*tau2[(i-1)]) / (sigma2[(i-1)] + n*tau2[(i-1)]))*Ybar), 
      Sigma=
        diag((sigma2[(i-1)]*tau2[(i-1)]) / (sigma2[(i-1)] + n * tau2[(i-1)]))
    )
    
    # MU
    mu[i] <- rnorm(n = 1, mean = (    
      (tau2[(i-1)] / (tau2[(i-1)] + k*sigma2u[(i-1)]))*mu0 + 
        ((k*sigma2u[(i-1)]) / (tau2[(i-1)] + k*sigma2u[(i-1)]))*mean(theta[(i-1),])), 
      sd = (
        sqrt((sigma2u[(i-1)]*tau2[(i-1)]) / (tau2[(i-1)] + k*sigma2u[(i-1)]))
      )
    )
    
    # SIGMA2
    sigma2[i] <- rinvgamma(n= 1, (ntot/2 + a[1]), ((1/2)*sum((data$Y- theta[i,])**2)) +b[1])
    
    # TAU2
    tau2[i] <- rinvgamma(n= 1, (k/2 + a[2]), (1/2)*sum((theta[i,] - mu[i])**2) + b[2])
    
    # SIMGA2u
    sigma2u[i] <- rinvgamma(n= 1, (1/2 + a[3]), (1/2)*(mu[i]-mu0)**2 + b[3])  
  }
  return(list(theta =theta[-1:-burn,], mu =mu[-1:-burn], 
              sigma2 =sigma2[-1:-burn], tau2 = tau2[-1:-burn],sigma2u = sigma2u[-1:-burn]))
  
}

# The simulation
sims <- hGibbs(burn=1000, niter=10000, data= data, a= a, b= b, mu0= mu0)


### Plottning
# Thetas.. 
results <- X[1000:10000, ]
theta <- results[,1:6]
par(mfrow= c(3,3))
for(i in 4:6){
  title <- paste('Theta', i)
  plot(theta[,i], type='l', main = title, ylab='',col=4,las=1)
  hist(theta[,i], main = title, prob = TRUE,las=1)
  lines(density(theta[,i]),col=4,las=1)
  plot(cumsum(theta[,i]) / 1:length(theta[,i]), type = 'l', main = title, ylab = '',las=1,col=4)
  
}

### Mu och sigmas
titles <- c('Mu','Sigma Square','Tau Square','Sigma_mu Square')
par(mfrow= c(1,3))
for(i in 8:8){
  title <- titles[i-6]
  plot(results[,i], type='l', main = title, ylab='',col=4,las=1)
  hist(results[,i], main = title, prob = TRUE,las=1)
  lines(density(results[,i]),col=4,las=1)
  plot(cumsum(results[,i]) / 1:length(results[,i]), type = 'l', main = title, ylab = '',las=1,col=4)
  
}



colnames(theta) <- c('Theta1', 'Theta2', 'Theta3', 'Theta4', 'Theta5', 'Theta6')
# A table of results
results <- cbind(theta, as.data.frame(sims[-1]))
Mean <- apply(results, 2, mean)
Sd <- apply(results, 2, sd)
CI <- HPDinterval(mcmc(results))

options(scipen=10)

t <-round(cbind(Mean, Sd, CI), digits=3)

xtable(t)

# Comparison with the MCMCpack --------------------------------------------


model <- MCMChregress(fixed= Y~ 1, random= ~1, group = 'Group', 
                      data = data, burnin = 1000, mcmc=10000, thin = 1, 
                      seed= NA, beta.start = 1500, sigma2.start = 1, 
                      Vb.start = 1, mubeta = 0, Vbeta = 1, 
                      r=1, R= .1, nu = 0.001, delta = 0.001)



#== MCMC analysis

# Graphics
pdf("Posteriors-MCMChregress.pdf")
plot(model$mcmc)
dev.off()

# Summary
summary(model$mcmc)

model$mcmc

