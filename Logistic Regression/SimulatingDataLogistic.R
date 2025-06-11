set.seed(1)

simulate_logistic_data <- function(n) {
  # Generate predictor variables
  X1 <- rnorm(n, mean = 0, sd = 1)  # Standard normal distribution
  X2 <- rnorm(n, mean = 0, sd = 1)
  
  # Define coefficients (intercept and slopes)
  beta_0 <- -0.5  # Intercept
  beta_1 <- 1.2   # Coefficient for X1
  beta_2 <- -0.8  # Coefficient for X2
  
  # Compute linear predictor
  linear_pred <- beta_0 + beta_1 * X1 + beta_2 * X2
  
  # Compute probability using logistic function
  prob <- 1 / (1 + exp(-linear_pred))
  
  # Generate binary response variable Y from Bernoulli distribution
  Y <- rbinom(n, size = 1, prob = prob)
  
  # Create a data frame
  data <- data.frame(X1 = X1, X2 = X2, Y = Y)
  
  return(data)
}

# Example: Generate dataset with 100 samples
simulated_data <- simulate_logistic_data(1e5)

X <- as.matrix(simulated_data[, c("X1", "X2")])
Y <- simulated_data$Y
Dx <- cbind(1, X)

save(X, Dx, Y, file = "SimulatedData.RData")