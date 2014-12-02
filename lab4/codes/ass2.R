#Data should be stored in y row by row so that an experiment
#with K factors and J reps will be stored in a KxJ matrix.
#The default choices for the hyperparameters should be changed to
#improve the performance of the algorithm

gibbs=function(y, istart = 50, ikeep = 200, mu0 = 0, s0 = 1, a1 = 0, b1 = 1, a2 = 0, 
	b2 = 1)
{
	K <- dim(y)[1]
	J <- dim(y)[2]
	onej <- rep(1, J)
	ybar <- apply(y, 1, mean)
	## 'apply' returns a vector or array or list of values obtained by 
	## applying a function to margins of an array or matrix
	## parameter 1 means row , 0 means column
	mu <- rep(0, (istart + ikeep))
	se <- rep(0, (istart + ikeep))
	st <- rep(0, (istart + ikeep))
	theta <- matrix(0, nrow=K , ncol = (istart + ikeep))
	mu[1] <- mu0
	se[1] <- b2/(2 * max(a2, 1))
	st[1] <- b1/(2 * max(a1, 1))
	theta[, 1] <- ybar
	for(i in 2:(istart + ikeep)) {
		stm <- st[i - 1]
		sem <- se[i - 1]
		mum <- mu[i - 1]
		thetam <- theta[, (i - 1)]
		mu[i] <- rnorm(1, (stm * mu0 + s0 * sum(thetam))/(stm + K * s0), sqrt((stm * s0)/(stm + K * s0)))
		st[i] <- (b1 + 0.5 * sum((thetam - mum) * (thetam - mum)))/ rgamma(1, a1 + K/2)
		se[i] <- (b2 + 0.5 * sum((y - thetam %o% onej)^2))/rgamma(1, a2 + (K * J)/2)
		theta1 <- rnorm(K, (sem * mum)/(J * stm + sem), sqrt((stm * sem)/(J * stm + sem)))
		theta[, i] <- theta1 + ((J * stm)/(J * stm + sem)) * ybar}
#We output the burn-in trials as well (for graphing purposes)
                out <- list(Theta = theta, Mu = mu, VarianceE = se, VarianceTheta = st)
	out
}
