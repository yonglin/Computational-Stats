library(SuppDists)
data_set <- c(1545,1440,1440,1520,1580,1540,1555,1490,1560,1495,
              1595,1550,1605,1510,1560,1445,1440,1595,1465,1545,
              1595,1630,1515,1635,1625,1520,1455,1450,1480,1445)

Y <- matrix(data_set,6,5)

#### initialize constants and parameters ####

N <- 20000
## lenght of chain
X <- matrix(0,N,10)
## the chain, a bivariate sample

Y_bar_1 <- mean(Y[1,])
Y_bar_2 <- mean(Y[2,])
Y_bar_3 <- mean(Y[3,])
Y_bar_4 <- mean(Y[4,])
Y_bar_5 <- mean(Y[5,])
Y_bar_6 <- mean(Y[6,])

n_i <- dim(Y)[2]
k <- dim(Y)[1]
n <- k * n_i

a_i <- 0.0001
b_i <- 0.0001

X[1,] <- c(Y_bar_1,Y_bar_2,Y_bar_3,Y_bar_4,Y_bar_5,Y_bar_6,1,1,1,1)

sum_sq = function(vec, c)
{
  temp <- (vec - c)^2
  sum_sq <- sum(temp)
  return (sum_sq)
}

sum_sq_matrix = function(M,c)
{
  temp <- 0
  for (i in 1:6)
  {
    temp <- temp + sum_sq(M[i,],c[i])
  }
  
  return (temp)
}


#### generate the chain ####
for (i in 2:N){
  P <- X[i-1, ]
  m <- P[7]*P[8] / (P[8] + n_i*P[9]) + P[9]*n_i*Y_bar_1 / (P[8] + n_i*P[9])
  s <- sqrt(P[8]*P[9] / (P[8] + n_i*P[9]))
  X[i,1] <- rnorm(1, m, s)
  X[i,2:10] <- X[i-1,2:10]
  
  P <- X[i, ]
  m <- P[7]*P[8] / (P[8] + n_i*P[9]) + P[9]*n_i*Y_bar_2 / (P[8] + n_i*P[9])
  s <- sqrt(P[8]*P[9] / (P[8] + n_i*P[9]))
  X[i, 2] <- rnorm(1, m, s)
  X[i,3:10] <- X[i-1,3:10]
  
  P <- X[i, ]
  m <- P[7]*P[8] / (P[8] + n_i*P[9]) + P[9]*n_i*Y_bar_3 / (P[8] + n_i*P[9])
  s <- sqrt(P[8]*P[9] / (P[8] + n_i*P[9]))
  X[i, 3] <- rnorm(1, m, s)
  X[i,4:10] <- X[i-1,4:10]
  
  P <- X[i, ]
  m <- P[7]*P[8] / (P[8] + n_i*P[9]) + P[9]*n_i*Y_bar_4 / (P[8] + n_i*P[9])
  s <- sqrt(P[8]*P[9] / (P[8] + n_i*P[9]))
  X[i, 4] <- rnorm(1, m, s)
  X[i,5:10] <- X[i-1,5:10]
  
  P <- X[i, ]
  m <- P[7]*P[8] / (P[8] + n_i*P[9]) + P[9]*n_i*Y_bar_5 / (P[8] + n_i*P[9])
  s <- sqrt(P[8]*P[9] / (P[8] + n_i*P[9]))
  X[i, 5] <- rnorm(1, m, s)
  X[i,6:10] <- X[i-1,6:10]
  
  P <- X[i, ]
  m <- P[7]*P[8] / (P[8] + n_i*P[9]) + P[9]*n_i*Y_bar_6 / (P[8] + n_i*P[9])
  s <- sqrt(P[8]*P[9] / (P[8] + n_i*P[9]))
  X[i, 6] <- rnorm(1, m, s)
  X[i,7:10] <- X[i-1,7:10]
  
  P <- X[i, ]
  m <- (mean(P[1:6]))*k*P[10] / (P[9] + k*P[10])
  s <- sqrt(P[9]*P[10] / (P[9] + k*P[10]))
  X[i, 7] <- rnorm(1, m, s)
  X[i,8:10] <- X[i-1,8:10]
  
  P <- X[i, ]
  m <- 0.5*n + a_i
  s <- (0.5*sum_sq_matrix(Y, P[1:6])) + b_i
  X[i, 8] <- rinvgamma(1, m, s)
  X[i,9:10] <- X[i-1,9:10]
  
  P <- X[i, ]
  m <- 0.5*k + a_i
  s <- (0.5* n_i *(sum_sq(P[1:6],P[7]))) + b_i
  X[i, 9] <- rinvgamma(1, m, s)
  X[i,10] <- X[i-1,10]
  
  P <- X[i, ]
  m <- 0.5 + a_i
  s <- 0.5*P[7]^2 + b_i
  X[i, 10] <- rinvgamma(1, m, s)
  
  
}

x <- X


hpd<-function(x,p){ 
  #generate an hpd set of level p, based 
  #on a sample x from the posterior 
  dx<-density(x) 
  md<-dx$x[dx$y==max(dx$y)] 
  px<-dx$y/sum(dx$y) 
  pxs<--sort(-px) 
  ct<-min(pxs[cumsum(pxs)< p]) 
  list(hpdr=range(dx$x[px>=ct]),mode=md) 
} 

dyes <- matrix(scan(sep=","),byrow=T,ncol=5)
1545, 1440, 1440, 1520, 1580
1540, 1555, 1490, 1560, 1495
1595, 1550, 1605, 1510, 1560
1445, 1440, 1595, 1465, 1545
1595, 1630, 1515, 1635, 1625
1520, 1455, 1450, 1480, 1445

dyes <- data.frame(cbind(as.factor(1:6),dyes))
dimnames(dyes)[[2]] <- c("Batch",paste("Y",1:5,sep=""))
dyes

dyes <- scan(sep=",")
1545, 1440, 1440, 1520, 1580
1540, 1555, 1490, 1560, 1495
1595, 1550, 1605, 1510, 1560
1445, 1440, 1595, 1465, 1545
1595, 1630, 1515, 1635, 1625
1520, 1455, 1450, 1480, 1445

dyes <- data.frame(Batch=rep(1:6,rep(5,6)),Y=dyes)
dyes$Batch <- as.factor(dyes$Batch)
dyes
