joint_proposal_SD_adapter_R <- function(accepted_variable, current_iteration, iteration_adapting_began,
                                        current_scaling_factor, mu_previous, current_parameter_values,
                                        current_covariance_matrix) {
  
  # Calculating the number of iterations since you started adapting, and from that the cooldown factor.
  # Basically this algorithm adapts the covariance matrix specifying the size and shape of your jumps.
  # Over time however, it changes the covariance matrix less and less, so that eventually you get an
  # unchanging covariance matrix. This is what the cooldown factor is.
  iterations_since_adapting_began <- current_iteration - iteration_adapting_began
  cooldown <- (iterations_since_adapting_began + 1)^-0.6
  
  # Calculating the various components
  new_correlation_matrix = ((1 - cooldown) * current_covariance_matrix) + (cooldown * ((t(current_parameter_values - mu_previous)) %*%  (current_parameter_values - mu_previous)))
  new_mu = ((1 - cooldown) * mu_previous) + (cooldown * current_parameter_values)
  log_new_scaling_factor = log(current_scaling_factor) + cooldown * (accepted_variable - 0.25)
  new_scaling_factor = exp(log_new_scaling_factor)
  new_covariance_matrix = new_scaling_factor * new_correlation_matrix
  
  # Outputting the relevant results
  output_list <- list()
  output_list[["New_Covariance_Matrix"]] = new_covariance_matrix
  output_list[["New_Mu"]] = new_mu
  output_list[["New_Scaling_Factor"]] = new_scaling_factor
  
  return(output_list)
  
  
}

adaptive_MCMC_running <- function(number_of_iterations, parameters_vector, sd_proposals, start_adaptation, end_adaptation) {
  
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
  
  # Components for Adaptive MCMC  
  Accepted <- 0
  current_scaling_factor <- 1
  current_mu <- matrix(nrow = 1, ncol = 5)
  current_mu[1, ] <- MCMC_output[1, ]
  current_covariance_matrix <- sigma
  
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
    
    proposed_parameter_values <- proposal_function(as.vector(MCMC_output[i, ]), current_covariance_matrix)
    proposed_posterior <- posterior_function(proposed_parameter_values, parameters_vector)
    
    likelihood_ratio = exp(proposed_posterior - current_posterior);
    
    if(runif(1) < likelihood_ratio) {
      
      MCMC_output[i + 1, ] <- proposed_parameter_values
      Acceptances[i] <- 1
      Logposterior_storage[i] <- proposed_posterior
      current_posterior <- proposed_posterior
      Accepted <- 1
      
    } else{
      
      MCMC_output[i + 1, ] <- MCMC_output[i, ]
      Acceptances[i] <- 0
      Logposterior_storage[i] <- current_posterior
      Accepted <- 0
      
    }
    
    Acceptance_Ratio[i] <- sum(Acceptances)/i
    
    if (i > start_adaptation & i < end_adaptation) {
      
      latest_parameter_values <- MCMC_output[i + 1, ]
      adapter_output <- joint_proposal_SD_adapter_R(Accepted, i, start_adaptation,
                                                   current_scaling_factor, current_mu, 
                                                   latest_parameter_values, current_covariance_matrix)
      
      current_mu <- adapter_output[["New_Mu"]];
      current_covariance_matrix <- adapter_output[["New_Covariance_Matrix"]];
      current_scaling_factor <- adapter_output[["New_Scaling_Factor"]];
      
    }
    
    if(i %% 1 == 0) {
      
      print(c("The iteration number is", i))
      print(c("The acceptance ratio is", Acceptance_Ratio[i]))
      print(current_mu)
      print(current_covariance_matrix)
      
    }
    
  }
  
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Posterior_Output"]] <- Logposterior_storage
  list[["Acceptance Ratio"]] <- Acceptance_Ratio
  return(list)
  
}


tic()
adaptive_MCMC_running(10, parameters_vector, sd_proposals, 2, 10)
toc()
    