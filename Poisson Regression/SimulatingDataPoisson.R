set.seed(1)

simulate_poisson_data <- function(n) {
  # Generate predictor variables
  X1 <- rnorm(n, mean = 0, sd = 2)
  X2 <- rnorm(n, mean = 0, sd = 2)
  
  # Define coefficients (intercept and slopes)
  beta_0 <- -0.5  # Intercept
  beta_1 <- 1.2   # Coefficient for X1
  beta_2 <- -0.8
  
  # Compute linear predictor
  Y_vec <- exp(beta_0 + beta_1 * X1 + beta_2 * X2)
  Y <- numeric(length(X1))
  for (i in 1:length(Y)) {
    Y[i] <- rpois(1, Y_vec[i])
  }
  
  # Create a data frame
  data <- data.frame(X1 = X1, X2 = X2, Y = Y)
  
  return(data)
}


# Example: Generate dataset with 100 samples
simulated_data <- simulate_poisson_data(1e5)

X <- as.matrix(simulated_data[, c("X1", "X2")])
Y <- simulated_data$Y
Dx <- cbind(1, X)

save(X, Dx, Y, file = "SimulatedData.RData")
glm(formula = Y ~ X, family = poisson)