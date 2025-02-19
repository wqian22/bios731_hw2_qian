# Define the negative log-likelihood function
neg_log_likelihood <- function(beta, X, Y) {
  eta <- beta[1] + beta[2] * X
  log_lik <- sum(Y * eta - log(1 + exp(eta))) 
  return(-log_lik)
}

# Define Newton's method
newton_method <- function(X, Y, tol = 1e-6, max_iter = 100) {
  beta <- c(0, 0)  # initial guess
  iter <- 0
  start_time <- Sys.time()
  
  while (iter < max_iter) {
    eta <- beta[1] + beta[2] * X
    p <- 1 / (1 + exp(-eta))
    
    W <- diag(p * (1 - p), n, n)
    X_mat <- cbind(1, X)
    
    grad <- t(X_mat) %*% (Y - p)
    Hessian <- -t(X_mat) %*% W %*% X_mat
    
    beta_new <- beta - solve(Hessian) %*% grad
    
    if (max(abs(beta_new - beta)) < tol) break
    beta <- beta_new
    iter <- iter + 1
  }
  
  end_time <- Sys.time()
  list(beta = as.vector(beta), time = end_time - start_time, iter = iter)
}


mm_algorithm <- function(X, Y, tol = 1e-6, max_iter = 100) {
  n <- length(Y)  # number of observations
  X_mat <- cbind(1, X)
  p <- ncol(X_mat)    # number of predictors (including intercept)
  
  beta <- c(0, 0)
  iter <- 0 
  start_time <- Sys.time()  # start timing
  
  while (iter < max_iter) {
    eta <- X_mat %*% beta 
    p_i <- 1 / (1 + exp(-eta))
    
    # Compute weight term based on minorization function
    weight <- exp(eta) / (1 + exp(eta))
    
    # Update each beta_j using the formula in 1.2(C)
    beta_new <- beta
    for (j in 1:p) {
      numerator <- sum(Y * X_mat[, j])
      denominator <- sum(weight * X_mat[, j] * exp(-p * X_mat[, j] * beta[j])) * exp(p * X_mat[, j] * beta[j])
      beta_new[j] <- log(numerator / denominator) / p
    }
    
    # Check for convergence
    if (max(abs(beta_new - beta)) < tol) break
    beta <- beta_new
    iter <- iter + 1
  }
  
  end_time <- Sys.time()  # end timing
  return(list(beta = beta, time = end_time - start_time, iter = iter))
}
