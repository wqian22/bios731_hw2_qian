library(here)
library(tidyverse)

###############################################################
## define or source functions used in code below
###############################################################

source(here("source", "01_simulate_data.R"))
source(here("source", "02_apply_methods.R"))

###############################################################
## set simulation design elements
###############################################################

nsim = 1
n = 200
beta0 = 1
beta1 = 0.3

###############################################################
## start simulation code
###############################################################

# generate a random seed for each simulated dataset
  seed <- floor(runif(nsim, 1, 900))
  results = vector("list", length = nsim)
  
  # use parallel computing for the simulations
  for(i in 1:nsim){

    set.seed(seed[i])
    
    ####################
    # simulate data
    simdata <- get_simdata(n = n, beta0 = beta0, beta1 = beta1)
    
    ####################
    # apply method
    newton_res <- newton_method(simdata$X, simdata$Y)
    mm_res <- mm_algorithm(simdata$X, simdata$Y)
    
    start_time <- Sys.time()
    glm_res <- glm(Y ~ X, data = simdata, family = binomial)
    end_time <- Sys.time()
    glm_time <- end_time - start_time
    glm_coefs <- coef(glm_res)
    glm_se <- sqrt(diag(vcov(glm_res)))
    
    start_time <- Sys.time()
    optim_res <- optim(par = c(0,0), fn = neg_log_likelihood, X = simdata$X, Y = simdata$Y, 
                       method = "BFGS", hessian = T)
    end_time <- Sys.time()
    optim_time <- end_time - start_time
    optim_coefs <- optim_res$par
    optim_se <- sqrt(diag(solve(optim_res$hessian)))
    
    ####################
    # Confidence intervals
    newton_ci <- cbind(newton_res$beta - 1.96 * glm_se, newton_res$beta + 1.96 * glm_se)
    mm_ci <- cbind(mm_res$beta - 1.96 * glm_se, mm_res$beta + 1.96 * glm_se)
    glm_ci <- cbind(glm_coefs - 1.96 * glm_se, glm_coefs + 1.96 * glm_se)
    optim_ci <- cbind(optim_coefs - 1.96 * optim_se, optim_coefs + 1.96 * optim_se)
    
    # Combine results into a table
    res <- data.frame(
      Method = c("Newton", "MM", "GLM", "BFGS"),
      Beta_0 = round(c(newton_res$beta[1], mm_res$beta[1], glm_coefs[1], optim_coefs[1]), 3),
      Beta_1 = round(c(newton_res$beta[2], mm_res$beta[2], glm_coefs[2], optim_coefs[2]), 3),
      CI_Beta_0 = paste0("(", round(newton_ci[1, 1], 3), ", ", round(newton_ci[1, 2], 3), ")"),
      CI_Beta_1 = paste0("(", round(newton_ci[2, 1], 3), ", ", round(newton_ci[2, 2], 3), ")"),
      Time = round(as.numeric(c(newton_res$time, mm_res$time, glm_time, optim_time)), 3),
      Iterations = c(newton_res$iter, mm_res$iter, glm_res$iter, optim_res$counts[1]/2)
    )
    
    results[[i]] <- res
  }
  
  results_df <- bind_rows(results)
  
  ####################
  # save results
  filename <- "results.RDA"
  save(results_df, file = here("results", filename))


