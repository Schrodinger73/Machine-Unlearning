library(mvtnorm)
library(rbenchmark)
library(ks)
library(FNN)
library(kableExtra)
library(dplyr)
library(coda)
library(DescTools)
library(MCMCglmm)
set.seed(1)

load("SimulatedData.RData")

g <- glm(Y ~ X, family = binomial)
g$coefficients

I <- which(Y == 1)
#Ie <- sample(I, size = 5e1, replace = FALSE)
Ie <- sample(I, size = 1e2, replace = FALSE)
Xe <- X[Ie,]
Ye <- Y[Ie]

Xr <- X[-Ie,]
Yr <- Y[-Ie]
gr <- glm(Yr ~ Xr, family = binomial)
gr$coefficients


log.f.real <- function(b) {
  s <- 0
  for (j in 1:length(Y)) {
    s <- s + as.numeric(Y[j]*t(Dx[j,]) %*% b - log(1 + exp(t(Dx[j,]) %*% b)))
  }
  return(s)
}

varMulti <- function(m) {
  mu <- rowMeans(m)
  V <- matrix(0, length(mu), length(mu))
  for (i in 1:length(m[1, ])) {
    V <- V + m[, i] %*% t(m[, i])
  }
  return(V/length(m[1, ]))
}

drawSamples <- function(n, h = 1) {
  M <- matrix(0, 3, n)
  a <- 0
  M[,1] <- rmvnorm(1, c(0, 0, 0), sigma = h*diag(1, 3, 3))
  for (i in 2:n) {
    y <- M[,(i - 1)] + rmvnorm(1, c(0, 0, 0), sigma = h*diag(1, 3, 3))
    log.alpha <- log.f.real(as.vector(y)) - log.f.real(as.vector(M[,(i - 1)]))
    if (log(runif(1)) < log.alpha) {
      a <- a + 1
      M[,i] <- y
    }
    else {
      M[,i] <- M[,(i - 1)]
    }
  }
  print(a/n)
  return(M)
}

# beta <- drawSamples(1e3, h = 5e-2)
# 
# for (i in 1:n) {
#   log.F.samples[[toString(beta[,i])]] <- 0
#   for (j in 1:length(Y)) {
#     log.F.samples[[toString(beta[,i])]] <- log.F.samples[[toString(beta[,i])]] + as.numeric(Y[j]*t(Dx[j,]) %*% beta[,i] - log(1 + exp(t(Dx[j,]) %*% beta[,i])))
#   }
# }
# 
# for (i in 1:n) {
#   H.samples[[toString(beta[,i])]] <- matrix(0, length(Dx[1,]), length(Dx[1,]))
#   for (j in 1:length(Y)) {
#     H.samples[[toString(beta[,i])]] <- H.samples[[toString(beta[,i])]] + as.numeric(exp(t(Dx[j,]) %*% beta[,i])/((1 + exp(t(Dx[j,])%*%beta[,i]))^2)) * (Dx[j,] %*% t(Dx[j,]))
#   }
# }

load("Model Approx to the Fourth Normal Random.Rdata")

approx <- function(b, k = 5) {
  D <- list()
  for (i in 1:length(beta[1,])) {
    D[[toString(beta[,i])]] <- norm(beta[,i] - b, type = "2")
  }
  sorted_keys <- names(sort(unlist(D)))[1:k]
  sorted_values <- sort(unlist(D))[1:k]
  s <- 0
  h <- matrix(0, length(Dx[1,]), length(Dx[1,]))
  d <- 0
  for (i in 1:k) {
    s <- s + log.F.samples[[sorted_keys[i]]]/sorted_values[i]^2
    h <- h + H.samples[[sorted_keys[i]]]/sorted_values[i]^2
    d <- d + 1/sorted_values[i]^2
  }
  L <- list()
  L[["log.f"]] <- s/d
  L[["H"]] <- h/d
  return(L)
}



proposal <- function(X, Y) {
  X_new <- cbind(1, X)
  beta.hat <- glm(Y ~ X, family = binomial)$coefficients
  H <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  De <- numeric(length(colnames(X_new)))
  for (i in 1:length(X_new[, 1])) {
    H <- H - exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat)))^2 * (X_new[i, ]) %*% t(X_new[i, ])
    De <- De + (exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat))) - Y[i]) * X_new[i, ]
  }
  L <- list()
  L[["H"]] <- H
  L[["D"]] <- De
  return(L)
}

proposal_new_approx <- function(X, Y, k = 5) {
  X_new <- cbind(1, X)
  beta.hat <- glm(Y ~ X, family = binomial)$coefficients
  A <- approx(beta.hat, k = k)
  H <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  H_real <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  De <- numeric(length(colnames(X_new)))
  for (i in 1:length(X_new[, 1])) {
    H <- H - exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat)))^2 * (X_new[i, ]) %*% t(X_new[i, ])
    De <- De + (exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat))) - Y[i]) * X_new[i, ]
  }
  #  for (i in 1:length(X_r_new[,1])) {
  #    H_real <- H_real + exp(as.numeric((X_r_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_r_new[i, ]) %*% beta.hat)))^2 * (X_r_new[i, ]) %*% t(X_r_new[i, ])
  #  }
  H_real <- H + A[["H"]]
  L <- list()
  L[["H Approx"]] <- H
  L[["H"]] <- H_real
  L[["D"]] <- De
  return(L)
}

proposal_mala_approx <- function(X, Y, k = 5) {
  X_new <- cbind(1, X)
  beta.hat <- glm(Y ~ X, family = binomial)$coefficients
  A <- approx(beta.hat, k = k)
  De <- numeric(length(colnames(X_new)))
  for (i in 1:length(X_new[, 1])) {
    De <- De + (exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat))) - Y[i]) * X_new[i, ]
  }
  #  for (i in 1:length(X_r_new[,1])) {
  #    H_real <- H_real + exp(as.numeric((X_r_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_r_new[i, ]) %*% beta.hat)))^2 * (X_r_new[i, ]) %*% t(X_r_new[i, ])
  #  }
  L <- list()
  L[["D"]] <- De
  return(L)
}

proposal_mala_approx_beta <- function(X, Y, B, k = 5) {
  X_new <- cbind(1, X)
  beta.hat <- B
  A <- approx(beta.hat, k = k)
  De <- numeric(length(colnames(X_new)))
  for (i in 1:length(X_new[, 1])) {
    De <- De + (exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat))) - Y[i]) * X_new[i, ]
  }
  #  for (i in 1:length(X_r_new[,1])) {
  #    H_real <- H_real + exp(as.numeric((X_r_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_r_new[i, ]) %*% beta.hat)))^2 * (X_r_new[i, ]) %*% t(X_r_new[i, ])
  #  }
  L <- list()
  L[["D"]] <- De
  return(L)
}

proposal_new_approx_beta <- function(X, Y, B, k = 5) {
  X_new <- cbind(1, X)
  beta.hat <- B
  A <- approx(beta.hat, k = k)
  H <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  H_real <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  De <- numeric(length(colnames(X_new)))
  for (i in 1:length(X_new[, 1])) {
    H <- H - exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat)))^2 * (X_new[i, ]) %*% t(X_new[i, ])
    De <- De + (exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat))) - Y[i]) * X_new[i, ]
  }
  #  for (i in 1:length(X_r_new[,1])) {
  #    H_real <- H_real + exp(as.numeric((X_r_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_r_new[i, ]) %*% beta.hat)))^2 * (X_r_new[i, ]) %*% t(X_r_new[i, ])
  #  }
  H_real <- H + A[["H"]]
  L <- list()
  L[["H Approx"]] <- H
  L[["H"]] <- H_real
  L[["D"]] <- De
  return(L)
}

proposal_new_real <- function(X, Y, Xr, Yr) {
  X_new <- cbind(1, X)
  X_r_new <- cbind(1, Xr)
  beta.hat <- glm(Y ~ X, family = binomial)$coefficients
  H <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  H_real <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  De <- numeric(length(colnames(X_new)))
  for (i in 1:length(X_new[, 1])) {
    H <- H - exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat)))^2 * (X_new[i, ]) %*% t(X_new[i, ])
    De <- De + (exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat))) - Y[i]) * X_new[i, ]
  }
  for (i in 1:length(X_r_new[,1])) {
    H_real <- H_real + exp(as.numeric((X_r_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_r_new[i, ]) %*% beta.hat)))^2 * (X_r_new[i, ]) %*% t(X_r_new[i, ])
  }
  L <- list()
  L[["H Approx"]] <- H
  L[["H"]] <- H_real
  L[["D"]] <- De
  return(L)
}

proposal_new_real_beta <- function(X, Y, Xr, Yr, B) {
  X_new <- cbind(1, X)
  X_r_new <- cbind(1, Xr)
  beta.hat <- B
  H <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  H_real <- matrix(0, nrow = length(colnames(X_new)), ncol = length(colnames(X_new)))
  De <- numeric(length(colnames(X_new)))
  for (i in 1:length(X_new[, 1])) {
    H <- H - exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat)))^2 * (X_new[i, ]) %*% t(X_new[i, ])
    De <- De + (exp(as.numeric((X_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_new[i, ]) %*% beta.hat))) - Y[i]) * X_new[i, ]
  }
  for (i in 1:length(X_r_new[,1])) {
    H_real <- H_real + exp(as.numeric((X_r_new[i,]) %*% beta.hat))/(1 + exp(as.numeric((X_r_new[i, ]) %*% beta.hat)))^2 * (X_r_new[i, ]) %*% t(X_r_new[i, ])
  }
  L <- list()
  L[["H Approx"]] <- H
  L[["H"]] <- H_real
  L[["D"]] <- De
  return(L)
}

Ha <- proposal_new_approx(Xe, Ye)
Hr <- proposal_new_real(X, Y, Xr, Yr)
solve(Ha$H)
solve(Hr$H)

log.f <- function(X, Y, beta) {
  s <- 0
  for (i in 1:length(Y)) {
    # print(i)
    # print(X[i,])
    # print(beta)
    # print(X[i,] %*% beta)
    # print(log(1 + exp(X[i, ] %*% beta)))
    s <- s + Y[i] * (X[i, ] %*% beta) - log(1 + exp(X[i, ] %*% beta))
  }
  #  print(s)
  return(s)
}

log.den.norm <- function(x, mean, sigma_inv) {
  n <- length(x)
  return(as.numeric(-n/2*log(2*pi) * log(det(sigma_inv))/2 - 1/2 * t(x - mean) %*% sigma_inv %*% (x - mean)))
} 

mcmc_nuevo_approx <- function(Xr, Yr, Xe, Ye, n = 1e3, h = 1e-6, k = 5, burnin = 0) {
  XE <- cbind(1, Xe)
  XR <- cbind(1, Xr)
  X <- rbind(Xe, Xr)
  Y <- append(Ye, Yr)
  XB <- cbind(1, X)
  beta.not <- glm(Y ~ X, family = binomial)$coefficients
  M <- matrix(0, length(XE[1,]), n + burnin)
  M[,1] <- beta.not
  #  prop <- proposal(Xe, Ye)
  prop.not <- proposal_new_approx(Xe, Ye, k = 10)
  H.not <- prop.not[["H"]]
  D.not <- prop.not[["D"]]
  print(H.not)
  print(D.not)
  H_inv.not <- solve(H.not)
  a <- 0
  for (i in 2:(n + burnin)) {
    print(i)
    beta.not <- as.vector(M[,(i - 1)])
    # prop <- proposal_new_approx_beta(Xe, Ye, beta.not, k = 10)
    # H.not <- prop[["H"]]
    # D <- prop[["D"]]
    # H_inv <- solve(H)
    # print(H_inv)
    # print(D)
    # print(beta.not + h*H_inv %*% D)
    beta.one <- as.vector(rmvnorm(1, mean = beta.not + h*H_inv.not %*% D.not, sigma = h^2*H_inv.not))
#    beta.one <- as.vector(rmvnorm(1, mean = beta.not + H_inv.not %*% D.not, sigma = h*diag(x = 1, length(beta.not))))
    prop.one <- proposal_new_approx_beta(Xe, Ye, beta.one)
    H.one <- prop.one[["H"]]
    D.one <- prop.one[["D"]]
    # print(H.one)
    # print(D.one)
    H_inv.one <- solve(H.one)
    #   print(M)
    # print("1 - ")
    # print(log.f(XE, Ye, beta.one))
    # print("2 - ")
    # print(log.f(XE, Ye, beta.not))
    L.one <- approx(beta.one, k = k)
    L.not <- approx(beta.not, k = k)
#    print(L.one)
#    print(L.not)
#    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] + dmvnorm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma = h*diag(x = 1, length(beta.not)), log = TRUE) - dmvnorm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma = h*diag(x = 1, length(beta.not)), log = TRUE)
#    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] +  log.den.norm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma_inv = diag(x = 1, length(beta.not))/h) - log.den.norm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma_inv = diag(x = 1, length(beta.not))/h)
#    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] + dmvnorm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma = h*H_inv.not, log = TRUE) - dmvnorm(beta.not, mean = beta.one + H_inv %*% D, sigma = h*H_inv, log = TRUE)
    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] + log.den.norm(beta.one, mean = beta.not + h*H_inv.not %*% D.not, sigma_inv = H.not/h^2) - log.den.norm(beta.not, mean = beta.one + h*H_inv.one %*% D.one, sigma_inv = H.one/h^2)
#    print(log.alpha)
    if (is.infinite(log.alpha)) {
      M[,i] <- beta.not
      next
    }
    
    if (is.na(log.alpha)) {
      M[,i] <- beta.not
      next
    }
    if (log(runif(1)) < log.alpha) {
      M[,i] <- beta.one
      if (i > burnin) {
        a <- a + 1
      }
      prop.not <- prop.one
      H.not <- prop.one[["H"]]
      D.not <- prop.one[["D"]]
      H_inv.not <- H_inv.one
    }
    else{
      M[,i] <- beta.not
    }
  }
  L = list()
  L[["M"]] = M[,((burnin + 1):(n + burnin))]
  L[["A"]] = a/n
  return(L)
}

mcmc_mala <- function(Xr, Yr, Xe, Ye, n = 1e3, h = 1e-6, k = 5, burnin = 0) {
  XE <- cbind(1, Xe)
  XR <- cbind(1, Xr)
  X <- rbind(Xe, Xr)
  Y <- append(Ye, Yr)
  XB <- cbind(1, X)
  beta.not <- glm(Y ~ X, family = binomial)$coefficients
  M <- matrix(0, length(XE[1,]), n + burnin)
  M[,1] <- beta.not
  #  prop <- proposal(Xe, Ye)
  prop.not <- proposal_mala_approx(Xe, Ye, k = 10)
  D.not <- prop.not[["D"]]
  a <- 0
  for (i in 2:(n + burnin)) {
    print(i)
    beta.not <- as.vector(M[,(i - 1)])
    # prop <- proposal_new_approx_beta(Xe, Ye, beta.not, k = 10)
    # H.not <- prop[["H"]]
    # D <- prop[["D"]]
    # H_inv <- solve(H)
    # print(H_inv)
    # print(D)
    # print(beta.not + h*H_inv %*% D)
    beta.one <- as.vector(rmvnorm(1, mean = beta.not + h/2*D.not, sigma = h*diag(1, length(beta.not), length(beta.not))))
#    beta.one <- as.vector(rmvnorm(1, mean = beta.not + H_inv.not %*% D.not, sigma = h*diag(x = 1, length(beta.not))))
    prop.one <- proposal_mala_approx_beta(Xe, Ye, beta.one)
    D.one <- prop.one[["D"]]
#    print(D.one)
    #   print(M)
    # print("1 - ")
    # print(log.f(XE, Ye, beta.one))
    # print("2 - ")
    # print(log.f(XE, Ye, beta.not))
    L.one <- approx(beta.one, k = k)
    L.not <- approx(beta.not, k = k)
    #    print(L.one)
    #    print(L.not)
    #    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] + dmvnorm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma = h*diag(x = 1, length(beta.not)), log = TRUE) - dmvnorm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma = h*diag(x = 1, length(beta.not)), log = TRUE)
    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] +  log.den.norm(beta.one, mean = beta.not + h/2 * D.not, sigma_inv = diag(x = 1, length(beta.not))/h) - log.den.norm(beta.not, mean = beta.one + h/2 * D.one, sigma_inv = diag(x = 1, length(beta.not))/h)
    #    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] + dmvnorm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma = h*H_inv.not, log = TRUE) - dmvnorm(beta.not, mean = beta.one + H_inv %*% D, sigma = h*H_inv, log = TRUE)
    #    log.alpha <- -(log.f(XE, Ye, beta.one)) + L.one[["log.f"]] + (log.f(XE, Ye, beta.not)) - L.not[["log.f"]] + log.den.norm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma_inv = H.not/h) - log.den.norm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma_inv = H.one/h)
#    print(log.alpha)
    if (is.infinite(log.alpha)) {
      M[,i] <- beta.not
      next
    }
    
    if (is.na(log.alpha)) {
      M[,i] <- beta.not
      next
    }
    if (log(runif(1)) < log.alpha) {
      M[,i] <- beta.one
      if (i > burnin) {
        a <- a + 1
      }
      prop.not <- prop.one
      D.not <- prop.one[["D"]]
    }
    else{
      M[,i] <- beta.not
    }
  }
  L = list()
  L[["M"]] = M[,((burnin + 1):(n + burnin))]
  L[["A"]] = a/n
  return(L)
}




drawSamplesReal <- function(n, h = 1, burnin = 0) {
  M <- matrix(0, 3, n + burnin)
  a <- 0
  Dr <- cbind(1, Xr)
  M[,1] <- rmvnorm(1, c(0, 0, 0), sigma = h*diag(1, 3, 3))
  for (i in 2:(n + burnin)) {
    print(i)
    y <- M[,(i - 1)] + rmvnorm(1, c(0, 0, 0), sigma = h*diag(1, 3, 3))
    log.alpha <- log.f(Dr, Yr, as.vector(y)) - log.f(Dr, Yr, as.vector(M[,(i - 1)]))
    if (log(runif(1)) < log.alpha) {
      if (i > burnin) {
        a <- a + 1
      }
      M[,i] <- y
    }
    else {
      M[,i] <- M[,(i - 1)]
    }
  }
  L = list()
  L[["M"]] = M[,((burnin + 1):(n + burnin))]
  L[["A"]] = a/n
  return(L)
}

XE <- cbind(1, Xe)
XR <- cbind(1, Xr)
X <- rbind(Xe, Xr)
Y <- append(Ye, Yr)
XB <- cbind(1, X)
beta.not <- glm(Y ~ X, family = binomial)$coefficients
#  prop <- proposal(Xe, Ye)
prop <- proposal_new_real(Xe, Ye, Xr, Yr)
H <- prop[["H"]]
D <- prop[["D"]]
print(H)
print(D)
H_inv <- solve(H)
beta.not + H_inv %*% D

mcmc_nuevo <- function(Xr, Yr, Xe, Ye, n = 1e3, h = 1e-6) {
  XE <- cbind(1, Xe)
  XR <- cbind(1, Xr)
  X <- rbind(Xe, Xr)
  Y <- append(Ye, Yr)
  XB <- cbind(1, X)
  beta.not <- glm(Y ~ X, family = binomial)$coefficients
  M <- matrix(0, length(XE[1,]), n)
  M[,1] <- beta.not
  #  prop <- proposal(Xe, Ye)
  prop.not <- proposal_new_real(Xe, Ye, Xr, Yr)
  H.not <- prop[["H"]]
  D.not <- prop[["D"]]
  print(H.not)
  print(D.not)
  H_inv.not <- solve(H)
  a <- 0
  for (i in 2:n) {
    print(i)
    beta.not <- as.vector(M[,(i - 1)])
#    prop <- proposal_new_real_beta(Xe, Ye, Xr, Yr, beta.not)
#    H <- prop[["H"]]
#    D <- prop[["D"]]
#    print(H)
#    print(D)
#    H_inv <- solve(H)
#    print(beta.not + h*H_inv %*% D)
    beta.one <- as.vector(rmvnorm(1, mean = beta.not + h*H_inv.not %*% D.not, sigma = h*H_inv.not))
#    beta.one <- as.vector(rmvnorm(1, mean = beta.not + H_inv %*% D, sigma = h*diag(x = 1, length(beta.not))))
    prop.one <- proposal_new_real_beta(Xe, Ye, Xr, Yr, beta.one)
    H.one <- prop.one[["H"]]
    D.one <- prop.one[["D"]]
    # print(H.one)
    # print(D.one)
    H_inv.one <- solve(H.one)

    # print(beta.not)
    # print(beta.one)
    # #   print(M)
    # print("1 - ")
    # print(log.f(XE, Ye, beta.one))
    # print("2 - ")
    # print(log.f(XE, Ye, beta.not))
#    log.alpha <- -(log.f(XE, Ye, beta.one)) + log.f(XB, Y, beta.one) + (log.f(XE, Ye, beta.not)) - log.f(XB, Y, beta.not) + dmvnorm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma = h*diag(x = 1, length(beta.not)), log = TRUE) - dmvnorm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma = h*diag(x = 1, length(beta.not)), log = TRUE)
#    log.alpha <- -(log.f(XE, Ye, beta.one)) + log.f(XB, Y, beta.one) + (log.f(XE, Ye, beta.not)) - log.f(XB, Y, beta.not) + log.den.norm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma_inv = diag(x = 1, length(beta.not))/h) - log.den.norm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma_inv = diag(x = 1, length(beta.not))/h)
#    log.alpha <- -(log.f(XE, Ye, beta.one)) + log.f(XB, Y, beta.one) + (log.f(XE, Ye, beta.not)) - log.f(XB, Y, beta.not) + dmvnorm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma = h*H_inv.not, log = TRUE) - dmvnorm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma = h*H_inv.one, log = TRUE)
    log.alpha <- -(log.f(XE, Ye, beta.one)) + log.f(XB, Y, beta.one) + (log.f(XE, Ye, beta.not)) - log.f(XB, Y, beta.not) + log.den.norm(beta.one, mean = beta.not + H_inv.not %*% D.not, sigma_inv = H.not/h) - log.den.norm(beta.not, mean = beta.one + H_inv.one %*% D.one, sigma = H.one/h)
#    print(log.alpha)
    if (is.infinite(log.alpha)) {
      M[,i] <- beta.not
      next
    }
    
    if (is.na(log.alpha)) {
      M[,i] <- beta.not
      next
    }
    if (log(runif(1)) < log.alpha) {
      M[,i] <- beta.one
      a <- a + 1
      prop.not <- prop.one
      H.not <- prop.one[["H"]]
      D.not <- prop.one[["D"]]
      H_inv.not <- H_inv.one
    }
    else{
      M[,i] <- beta.not
    }
  }
  print(a/n)
  L = list()
  L[["M"]] = M
  L[["A"]] = a/n
  return(L)
}
# 2e-1
ta <- system.time(maF <- mcmc_nuevo_approx(Xr, Yr, Xe, Ye, n = 1000, h = 1e-1, k = 200, burnin = 200))
print("Approx")
# mnF <- mcmc_nuevo(Xr, Yr, Xe, Ye, n = 5e2, h = 5e-5)
# print("Unlearning")
# 6e-8
tm <- system.time(mmF <- mcmc_mala(Xr, Yr, Xe, Ye, n = 1000, h = 1e-7, k = 200, burnin = 200))
print("MALA")
# 0.00015
tr <- system.time(mrF <- drawSamplesReal(1000, h = 0.0001, burnin = 400))
print("Real")
ma <- maF$M
#mn <- mnF$M
mm <- mmF$M
mr <- mrF$M
rowMeans(ma)
#rowMeans(mn)
rowMeans(mm)
rowMeans(mr)
gr$coefficients
g$coefficients

b <- benchmark(mcmc_nuevo_approx(Xr, Yr, Xe, Ye, n = 1e2, h = 5e-5, k = 200), mcmc_mala(Xr, Yr, Xe, Ye, n = 1e2, k = 200, h = 8e-7), drawSamplesReal(1e2, h = 5e-2), replications = 1)
b
kable(b, "latex")

rowMeans(ma)
rowMeans(mm)
rowMeans(mr)
gr$coefficients
g$coefficients



k <- 150  # Number of neighbors
densities_a <- knn.dist(t(ma), k = k) %>% rowMeans()
#densities_n <- knn.dist(t(mn), k = k) %>% rowMeans()
densities_m <- knn.dist(t(mm), k = k) %>% rowMeans()
densities_r <- knn.dist(t(mr), k = k) %>% rowMeans()

# Find mode (highest density point)
mode_index_a <- which.min(densities_a)
mode_point_a <- t(ma)[mode_index_a, ]

# mode_index_n <- which.min(densities_n)
# mode_point_n <- t(mn)[mode_index_n, ]

mode_index_m <- which.min(densities_m)
mode_point_m <- t(mm)[mode_index_m, ]

mode_index_r <- which.min(densities_r)
mode_point_r <- t(mr)[mode_index_r, ]

rowMeans(ma)
rowMeans(mm)
rowMeans(mr)

varMulti(ma)
varMulti(mm)
varMulti(mr)

posterior.mode(mca)
posterior.mode(mcm)
posterior.mode(mcr)
maF$A
mmF$A
mrF$A
gr$coefficients
g$coefficients


mca <- mcmc(t(ma))
mcm <- mcmc(t(mm))
mcr <- mcmc(t(mr))

geweke.diag(mca)
geweke.diag(mcm)
geweke.diag(mcr)

autocorr(mca)
autocorr(mcm)
autocorr(mcr)

plot(mca)
plot(mcm)
plot(mcr)

rowMeans(ma)
rowMeans(mm)
rowMeans(mr)

varMulti(ma)
varMulti(mm)
varMulti(mr)

posterior.mode(mca)
posterior.mode(mcm)
posterior.mode(mcr)

maF$A
mmF$A
mrF$A
gr$coefficients
g$coefficients

min(effectiveSize(mca))
min(effectiveSize(mcm))
min(effectiveSize(mcr))

min(effectiveSize(mca)/ta[3])
min(effectiveSize(mcm)/tm[3])
min(effectiveSize(mcr)/tr[3])

geweke.diag(mca)
geweke.diag(mcm)
geweke.diag(mcr)
