get_simdata = function(n, beta0, beta1){
  X <- rnorm(n, mean = 0, sd = 1)
  P_Y <- 1 / (1 + exp(-(beta0 + beta1 * X)))
  Y <- rbinom(n, size = 1, prob = P_Y)

  tibble(
    X = X,
    Y = Y
  )

}



