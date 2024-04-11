setwd("~/Desktop/S2/X4MS040")

library(coda)

x <- c(0, 10, 33, 100, 333, 1000)
i1 <- c(15, 16, 16, 27, 33, 20)
i2 <- c(21, 18, 26, 41, 38, 27)
i3 <- c(29, 21, 33, 69, 41, 42)
y <- structure(c(15, 16, 16, 27, 33, 20, 21, 18, 26, 41, 38, 27, 29, 
                 21, 33, 60, 41, 42), .Dim = c(6, 3))


GIBBS <- function(x, y, nchain, init_abg, init_tau, sd_abg_prop){
  chain_abg <- matrix(NA, nchain+1, 3) # chaine pour alpha beta gamma
  chain_abg[1,] <- init_abg
  chain_tau <- matrix(NA, nchain+1, 1) # chaine pour tau
  chain_tau[1] <- init_tau
  
  # hyperparamètres
  a <- 1/1000
  b <- 1/1000
  sigma <- 100000
  
  
  for (i in 1:nchain){
    # MAJ tau
    #tau_candidate <- rnorm(1,mean=chain_tau[i], sd=sd_tau_prop)
    # Lambda nécessaire à la construction des paramètres tau et MAJ mu
    lambda <- rnorm(3*6, mean=0, sd=1/chain_tau[i])
    lambda <- matrix(lambda, nrow=6, ncol=3)
    #lambda_candidate <- rnorm(3*6, mean=0, sd=1/tau_candidate)
    #lambda_candidate <- matrix(lambda_candidate, nrow=6, ncol=3)
    
    #top <- log(rgamma(1, a + (1/2)*3*6, b + (1/2)*sum(lambda_candidate**2)))
    #bottom <-log(rgamma(1, a + (1/2)*3*6, b + (1/2)*sum(lambda**2)))
    #ratio <- top - bottom
    #u <- log(runif(1))
    
   # if (u < ratio){
    #  chain_tau[i+1] <- tau_candidate
    #} else {
    #  chain_tau[i+1] <- chain_tau[i]
    #}
    
    chain_tau[i+1] <- rgamma(1, a + (1/2)*3*6, b + (1/2)*sum(lambda**2))
    
    # MAJ alpha beta gamma
    lambda <- rnorm(3*6, mean=0, sd=1/chain_tau[i+1])
    lambda <- matrix(lambda, nrow=6, ncol=3)
    # alpha
    alpha_candidate <- rnorm(1, mean = chain_abg[i,1], sd = sd_abg_prop[1])
    mu_candidate <- exp(alpha_candidate + chain_abg[i,2]*log(x + 10) +
                          chain_abg[i,3]*x + lambda)
    logmu_candidate <- alpha_candidate + chain_abg[i,2]*log(x + 10) +
      chain_abg[i,3]*x + lambda
    mu <- exp(chain_abg[i,1] + chain_abg[i,2]*log(x + 10) +
                chain_abg[i,3]*x + lambda)
    logmu <- chain_abg[i,1] + chain_abg[i,2]*log(x + 10) +
      chain_abg[i,3]*x + lambda
    top <- (alpha_candidate**2)/(2*sigma) + sum(mu_candidate) - sum(y*logmu_candidate)
    bottom <- (chain_abg[i,1]**2)/(2*sigma) + sum(mu) - sum(y*logmu)
    ratio <- top - bottom
    u <- log(runif(1))
    
    if (u < ratio){
      chain_abg[i+1,1] <- alpha_candidate
    } else {
      chain_abg[i+1,1] <- chain_abg[i,1]
    }
    
    # Beta 
    beta_candidate <- rnorm(1, mean = chain_abg[i,2], sd = sd_abg_prop[2])
    mu_candidate <- exp(chain_abg[i+1,1] + beta_candidate*log(x + 10) +
                          chain_abg[i,3]*x + lambda)
    logmu_candidate <- chain_abg[i+1,1] + beta_candidate*log(x + 10) +
      chain_abg[i,3]*x + lambda
    mu <- exp(chain_abg[i+1,1] + chain_abg[i,2]*log(x + 10) +
                chain_abg[i,3]*x + lambda)
    logmu <- chain_abg[i+1,1] + chain_abg[i,2]*log(x + 10) +
      chain_abg[i,3]*x + lambda
    top <- (beta_candidate**2)/(2*sigma) + sum(mu_candidate) - sum(y*logmu_candidate)
    bottom <- (chain_abg[i,2]**2)/(2*sigma) + sum(mu) - sum(y*logmu)
    ratio <- top - bottom
    u <- log(runif(1))
    
    if (u < ratio){
      chain_abg[i+1,2] <- beta_candidate
    } else {
      chain_abg[i+1,2] <- chain_abg[i,2]
    }
    
    
    # Gamma
    gamma_candidate <- rnorm(1, mean = chain_abg[i,3], sd = sd_abg_prop[3])
    mu_candidate <- exp(chain_abg[i+1,1] + chain_abg[i+1,2]*log(x + 10) +
                          gamma_candidate*x + lambda)
    log_mucandidate <- chain_abg[i+1,1] + chain_abg[i+1,2]*log(x + 10) +
      gamma_candidate*x + lambda
    logmu <- chain_abg[i+1,1] + chain_abg[i+1,2]*log(x + 10) +
      chain_abg[i,3]*x + lambda
    mu <- exp(chain_abg[i+1,1] + chain_abg[i+1,2]*log(x + 10) +
                chain_abg[i,3]*x + lambda)
    top <- (gamma_candidate**2)/(2*sigma) + sum(mu_candidate) - sum(y*logmu_candidate)
    bottom <- (chain_abg[i,3]**2)/(2*sigma) + sum(mu) - sum(y*logmu)
    ratio <- top - bottom
    u <- log(runif(1))
    
    if (u < ratio){
      chain_abg[i+1,3] <- gamma_candidate
    } else {
      chain_abg[i+1,3] <- chain_abg[i,3]
    }
    
  }
  
  return(list(tau = chain_tau, abg = chain_abg))
  
}

result <- GIBBS(x, y, 10000, c(0,0,0), 1, c(1,1,1))


plot(mcmc(result$tau[1000:10001]))
plot(mcmc(result$abg[1000:10001,1]))
plot(mcmc(result$abg[1000:10001,2]))
plot(mcmc(result$abg[1000:10001,3]))

