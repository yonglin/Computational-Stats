# The second derivative f2 of f
fexp <- expression(log(1+2*x) / (1+x^2) )
f1exp <- D(fexp,"x")
f2exp <- D(f1exp,"x")

f <- function(x){
  eval(fexp)
  }
f1 <- function(x){
  eval(f1exp)
  }
f2 <- function(x){
  eval(f2exp)
  }
# Initial values
x <- 0
itr <- 40
h <- c(1:itr)
netraj <- c(1:itr)
neftraj <- c(1:itr)
# Newton’s method
for(i in 1:itr){
  netraj[i] <- x
  neftraj[i] <- f(x)
  h[i] <- f1(x)/f2(x)
  x <- x - h[i]

}

## Plotting
curve(f,from=0,to=10,las=1,col=2)
lines(netraj,neftraj,col=4,lty=2)
points(x,f(x),col=4)
text(x,f(x),round(x,5),pos=1)