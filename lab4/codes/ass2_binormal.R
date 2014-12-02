#### initialize constants and parameters ####

N <- 5000
## lenght of chain
burn <- 1000
## burn-in length
X <- matrix(0,N,2)
## the chain, a bivariate sample

rho <- -0.75
## corelation
mu_1 <- 0
mu_2 <- 2
sigma_1 <- 1
sigma_2 <- 0.5

s_1 <- sqrt(1-rho^2)*sigma_1
s_2 <- sqrt(1-rho^2)*sigma_2

#### generate the chain ####

X[1,] <- c(mu_1, mu_2)

for (i in 2:N){
  x_2 <- X[i-1,2]
  m_1 <- mu_1 + rho * (x_2 - mu_2) * sigma_1/sigma_2
  X[i, 1] <- rnorm(1, m_1, s_1)
  
  x_1 <- X[i, 1]
  m_2 <- mu_2 + rho * (x_1 - mu_1) * sigma_2/sigma_1
  X[i,2] <- rnorm(1, m_2, s_2)
}

b <- burn +1
x <- X

plot(x, main="", cex=.5, xlab=bquote(X[1]),ylab=bquote(X[2]), ylim=range(x[,2]),las=1, col=4)
