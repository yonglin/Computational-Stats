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