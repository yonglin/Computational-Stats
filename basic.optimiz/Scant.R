# The second derivative f2 of f
fexp <- expression(log(1+2*x) / (1+x^2) )
f1exp <- D(fexp,"x")


f <- function(x){
  eval(fexp)
}
f1 <- function(x){
  eval(f1exp)
}

# Initial values
x1 <- 0.2
x <- 0.5
itr <- 40
##h <- c(1:itr)
netraj <- c(1:itr)
neftraj <- c(1:itr)
netraj[1] <- x1
netraj[2] <- x
neftraj[1] <- f(x1)
neftraj[2] <- f(x)
# Newton’s method
for(i in 3:itr){
  ## The relative convergence criterion threshold delta
  delta <- (netraj[i-1]-netraj[i-2])/netraj[i-2]
  if (delta>0.01){
    x <- x - (f1(netraj[i-1])*(netraj[i-1]-netraj[i-2]))/(f1(netraj[i-1])-f1(netraj[i-2]))
    netraj[i] <- x
    neftraj[i] <- f(x)
  }
  else {break}
  
}
f(x)
## Plotting
curve(f,from=0,to=10,las=1,col=2)
lines(netraj,neftraj,col=4,lty=2)
points(x,f(x),col=4)
text(x,f(x),round(x,5),pos=1)