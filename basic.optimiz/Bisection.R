## The derivative f1 of f
fexp <- expression(log(1+2*x) / (1+x^2) )
f1exp <- D(fexp,"x")
f <- function(x){
  eval(fexp)
  }
f1 <- function(x){
  eval(f1exp)
  }

## Initial values
a <- 0
b <- 6
x <- a+(b-a)/2
itr <- 40
bitraj <- c(1:itr)
biftraj <- c(1:itr)

## Bisection
for (i in 1:itr){
  bitraj[i] <- x
  biftraj[i] <- f(x)
  if (f1(a)*f1(x) < 0)
    {b = x}
  else 
    {a = x}
  x = a+(b-a)/2
}

## Plotting
curve(f,from=0,to=10,las=1,col=2)
lines(bitraj,biftraj,col=4,lty=2)
points(x,f(x),col=4)
text(x,f(x),round(x,5),pos=1)
  
