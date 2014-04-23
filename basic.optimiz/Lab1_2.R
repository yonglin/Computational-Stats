## Optimization of Himmelblau's function using Newton's method
## 
## Initial values
x <- c(0,0)
itr <- 40
x.values <- matrix(0,itr+1,2)
x.values[1,] <- x

x_new <- c(-2,-2)
itr <- 40
x_new.values <- matrix(0,itr+1,2)
x_new.values[1,] <- x_new

## Objective function and derivatives
f <- function(x){((((x[1]^2)+x[2]-11)^2)+(x[1]+(x[2]^2)-7)^2)}

f.prime <- function(x){
  f.prime.da <- ((4*x[1]^3)+(4*x[1]*x[2])-(42*x[1])+(2*x[2]^2)-14)
  f.prime.db <- ((2*x[1]^2)-(26*x[2])-22+(4*x[1]*x[2])+(4*x[2]^3))
  out <- matrix(c(f.prime.da,f.prime.db),ncol=1)
  return(out)
}
f.2prime <- function(x){
  f.2prime.da2 <- ((12*x[1]^2)+(4*x[2])-42)
  f.2prime.db2 <- ((12*x[2]^2)+(4*x[1])-26)
  f.2prime.dadb <- (4*(x[1]+x[2]))
  out <- matrix(c(f.2prime.da2,f.2prime.dadb,
                  f.2prime.dadb,f.2prime.db2),nrow=2, byrow=TRUE)
  return(out)
}
## Newton's method
for(i in 1:itr){
  x <- x - solve(f.2prime(x))%*%f.prime(x)
  x.values[i+1,] <- x
}
for(i in 1:itr){
  x_new <- x_new - solve(f.2prime(x_new))%*%f.prime(x_new)
  x_new.values[i+1,] <- x_new
}

## Output
x	
f(x) 		
f.prime(x) 

x_new  
f(x_new) 		
f.prime(x_new) 

## Plot 
z = matrix(0,101,101)
x1.max = max(4.5,ceiling(max(x.values[,1])))
x1.min = min(-2,floor(min(x.values[,1])))
x2.max = max(3,ceiling(max(x.values[,2])))
x2.min = min(-3,floor(min(x.values[,2])))
x1 = seq(x1.min,x1.max,length=101)
x2 = seq(x2.min,x2.max,length=101)

for(i in 1:101){
  for(j in 1:101){
    z[i,j] = f(c(x1[i],x2[j]))
  }
}
contour(x1,x2,z,nlevels=20,col='blue',las=1)

for(i in 1:itr){
  segments(x_new.values[i,1],x_new.values[i,2],x_new.values[i+1,1],
           x_new.values[i+1,2],lty=2,col=2)
}

for(i in 1:itr){
  segments(x.values[i,1],x.values[i,2],x.values[i+1,1],
           x.values[i+1,2],col=1)
}

points(x[1],x[2],pch=25)
points(x_new[1],x_new[2],col=2,pch=25)