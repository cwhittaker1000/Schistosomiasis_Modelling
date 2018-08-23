Estimate_Schisto_Data <- Ready_Schisto 
Estimate_Schisto_Data$epg_base <- round(Estimate_Schisto_Data$epg_base)

likelihood_function_epg <- function(parameters) {
  
  likelihood <- 0
  #vector_liks <- vector(length = 204)
  for (i in 1:length(Ready_Schisto[, 1])) {
    #vector_liks[i] <- dnbinom(Estimate_Schisto_Data$epg[i], size = parameters[1], mu = parameters[2], log = T)
    likelihood <- likelihood + dnbinom(Estimate_Schisto_Data$epg[i], size = parameters[1], mu = parameters[2], log = T)
  }
  
  return(-likelihood)
  
}

parameters <- c(0.4063612, 1012.0310537)
fitted <- optim(parameters, fn = likelihood_function_epg, method = "BFGS")
                #lower = c(1, 0.001, upper = c(50, 2))
our_negbinom <- rnbinom(100000, size = fitted$par[1], mu = fitted$par[2])
density <- density(our_negbinom)
plot(density$x, density$y, type = "l", xlim = c(0, 10000))


likelihood_function_IgE <- function(parameters) {
  
  likelihood <- 0

  for (i in 1:length(Ready_Schisto[, 1])) {
    
    likelihood <- likelihood + dbeta(Ready_Schisto$TAL1_IgE_base[i], shape1 = parameters[1], shape2 = parameters[2], log = T)
  
  }
  return(-likelihood)
}




likelihood_function_IgE(c(1, 2))

par <- c(1, 10)
fitted_IgE <- optim(par, fn = likelihood_function_IgE, method = "L-BFGS-B",
                    lower = c(0.05, 1), upper = c(5, 50))
fitted_IgE$par

?dbeta


beta<- rnbinom(100000, size = 10, mu = 10)
density_beta<-density(beta)
plot(density_beta$x, density_beta$y, axes=T, main="", type = "l")

hist(Ready_Schisto$TAL1_IgE_base, xlim=c(0,1), breaks=20)
par(new=T)
plot(density_beta$x, density_beta$y, xlim=c(0,1), axes=T, main="", type = "l")
