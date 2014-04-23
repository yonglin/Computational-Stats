## Function from Givens and Heoting

fexp <- expression(log(1+2*x) / (1+x^2) )

f <- function(x){
  eval(fexp)
  }

curve(f, from=0,to=10,col=2,las=1)

max <- optimize(f,c(0,6),maximum=T)

points(max$maximum,max$objective,col=4)

text(max$maximum,max$objective,round(max$maximum,5),pos=1)
