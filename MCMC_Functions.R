source("Model_Driver.R")
library(MASS)

# Prior Function ## double check fact it returns a named vector
prior_function <- function(parameters_vector) {
  
  R0 <- parameters_vector[1]
  worm_lifespan <- parameters_vector[6]
  immunity_gain <- parameters_vector[14]
  immunity_loss <- parameters_vector[15]
  egg_prod_rate <- parameters_vector[16]
  
  R0_prior_value <- dunif(R0, min = 0.5, max = 4.5, log = T)
  worm_lifespan_prior_value <- dunif(worm_lifespan, min = 0.05, max = 0.5, log = T)  # change this so that it isn't the rate but the actual lifespan
  immunity_gain_prior_value <- dunif(immunity_gain, min = 0, max = 1, log = T)
  immunity_loss_prior_value <- dunif(immunity_loss, min = 0, max = 1, log = T)
  egg_prod_rate_prior_value <- dunif(egg_prod_rate, min = 2, max = 8, log = T)
  
  prior_output <- R0_prior_value + worm_lifespan_prior_value + immunity_gain_prior_value + immunity_loss_prior_value + egg_prod_rate_prior_value
  return(prior_output)

}

prior_function(parameters)

proposed_parameter_values <- c(3, 0.1, 0.1, 0.1, 6)

#Likelihood Function
loglikelihood_function <- function(proposed_parameter_values, parameters_vector) {
  
  # Input proposed parameters into parameters vector
  parameters_vector["R0"] <- proposed_parameter_values[1]
  parameters_vector["mu.w"] <- proposed_parameter_values[2]
  parameters_vector["gain.imm"] <- proposed_parameter_values[3]
  parameters_vector["loss.imm"] <- proposed_parameter_values[4]
  parameters_vector["eggrate"] <- proposed_parameter_values[5]
  
  # Running the model
  model_output <- runmod(parameters_vector)
  processed_output <- extract(model_output, parameters_vector)
  
  # Extracting the relevant model output
  age_points <- c(8, 13, 18, 23, 28, 33, 38, 43, 48)
  epg_specific_ages <- processed_output$epg[age_points]
  IgE_specific_ages <- processed_output$IgE[age_points]
  
  # Get the mean of the data for each age group 
  mean_epg_age_group <- round(aggregate(epg_base ~ agecat, data=Ready_Schisto, FUN="mean")[, 2])
  mean_IgE_age_group <- aggregate(TAL1_IgE_base ~ agecat, data=Ready_Schisto, FUN="mean")[, 2]
  
  # Calculate the loglikelihoods 
  loglik_epg <- sum(dnbinom(mean_epg_age_group, size = 0.4064213, mu = epg_specific_ages, log = T))
  loglik_IgE <- sum(dbeta(mean_IgE_age_group, shape1 = 1, shape2 = 20, log = T))
  
  # Return overall loglikelihood
  overall_loglikelihood <- loglik_epg + loglik_IgE
  return(overall_loglikelihood)
}

parameters_vector <- c(R0=3,a1=0.9,a2=15, a3=0.6, a4=0.1, mu.w=0.2, mu.h=0.04,c=88.95, beta=0.05, d=15 , 
                       nw = 5, na=50, amax=50, gain.imm= 0.1, loss.imm = 0.1, eggrate=6)
proposed_parameter_values <- c(3, 0.1, 0.1, 0.1, 6)
loglik <- loglikelihood_function(proposed_parameter_values, parameters_vector)

# Proposal Function
proposal_function <- function(current_parameters, covariance_matrix) {
  
  proposed_parameter_values <- mvrnorm(1, current_parameters, covariance_matrix)
  
  for (i in 1:length(proposed_parameter_values)) {
    if (proposed_parameter_values[i] <= 0) {
      proposed_parameter_values[i] <- 0.01
    }
  }
  
  return(proposed_parameter_values)
  
}  

Sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = T)
proposal_function(c(5, 10), Sigma)

# Posterior Function
posterior_function <- function(proposed_parameter_values, parameters_vector){
  
  posterior <- loglikelihood_function(proposed_parameter_values, parameters_vector) + prior_function(parameters_vector)
  return(posterior)
  
}

posterior_function(proposed_parameter_values, parameters_vector)
prior_function(parameters_vector)
loglikelihood_function(proposed_parameter_values, parameters_vector)

# MCMC Function
MCMC_running <- function(number_of_iterations, parameters_vector, sd_proposals) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol = 5)
  MCMC_output[1, ] <- c(parameters_vector[1], parameters_vector[6], parameters_vector[14], parameters_vector[15], parameters_vector[16])
  colnames(MCMC_output) = c("R0", "mu.W", "gain.Imm", "loss.Imm", "EggRate")
  Acceptances <- vector(length = number_of_iterations + 1)
  Acceptance_Ratio <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  # Generating the Covariance Matrix
  sigma <- matrix(0, nrow = 5, ncol = 5)
  for (i in 1:5) {
    for (j in 1:5) {
      if (i == j) {
        sigma[i, j] <- sd_proposals[i]
      }
    }
  }
  
  # Calculating the Initial Posterior
  fake_proposed_values <- vector(length = 5)
  fake_proposed_values[1] <- parameters_vector[1]
  fake_proposed_values[2] <- parameters_vector[6]
  fake_proposed_values[3] <- parameters_vector[14]
  fake_proposed_values[4] <- parameters_vector[15]
  fake_proposed_values[5] <- parameters_vector[16]
  current_posterior <- posterior_function(fake_proposed_values, parameters_vector)
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    proposed_parameter_values <- proposal_function(MCMC_output[i, ], sigma)
    proposed_posterior <- posterior_function(proposed_parameter_values, parameters_vector)
    
    likelihood_ratio = exp(proposed_posterior - current_posterior);
    
    if(runif(1) < likelihood_ratio) {
      
      MCMC_output[i + 1, ] <- proposed_parameter_values
      Acceptances[i] <- 1
      Logposterior_storage[i] <- proposed_posterior
      current_posterior <- proposed_posterior
      
    } else{
      
      MCMC_output[i + 1, ] <- MCMC_output[i, ]
      Acceptances[i] <- 0
      Logposterior_storage[i] <- current_posterior
      
    }
    
    Acceptance_Ratio[i] <- sum(Acceptances)/i
    
    if(i %% 10 == 0) {
      
      print(c("The iteration number is", i))
      print(c("The acceptance ratio is", Acceptance_Ratio[i]))

    }
    
  }
  
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Posterior_Output"]] <- Logposterior_storage
  return(list)
  
}

library(tictoc)
sd_proposals <- c(0.02, 0.005, 0.003, 0.003, 0.2)
tic()
output <- MCMC_running(2, parameters_vector, sd_proposals)
toc()






