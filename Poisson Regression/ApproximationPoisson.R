library(mvtnorm)
n <- 9e2
#beta <- drawSamples(n)
log.F.samples <- list()
H.samples <- list()

load("SimulatedData.RData")


log.f.real <- function(b) {
  s <- 0
  for (j in 1:length(Y)) {
    s <- s + as.numeric(Y[j]*t(Dx[j,]) %*% b - exp(t(Dx[j,]) %*% b))
  }
  return(s)
}

drawSamples <- function(n, h = 1) {
  M <- matrix(0, 3, n)
  a <- 0
  M[,1] <- rmvnorm(1, c(0, 0, 0), sigma = h*diag(1, 3, 3))
  for (j in 1:length(Y)) {
    H.samples[[toString(M[,1])]] <- H.samples[[toString(as.vector(M[,1]))]] + as.numeric(exp(t(Dx[j,]) %*% as.vector(M[,1]))) * Dx[j,] %*% t(Dx[j,])
  }
  for (i in 2:n) {
    print(i)
    y <- M[,(i - 1)] + rmvnorm(1, c(0, 0, 0), sigma = h*diag(1, 3, 3))
    k <- log.f.real(as.vector(y))
    l <- log.f.real(as.vector(M[, (i - 1)]))
    log.alpha <- k - l
    log.F.samples[[toString(as.vector(y))]] <- k
    log.F.samples[[toString(as.vector(M[,(i - 1)]))]] <- l
    log.alpha <- k - l
    if (log(runif(1)) < log.alpha) {
      a <- a + 1
      M[,i] <- y
      for (j in 1:length(Y)) {
        H.samples[[toString(as.vector(y))]] <- H.samples[[toString(as.vector(y))]] + as.numeric(exp(t(Dx[j,]) %*% as.vector(y))) * Dx[j,] %*% t(Dx[j,])
      }
    }
    else {
      M[,i] <- M[,(i - 1)]
    }
  }
  print(a/n)
  return(M)
}
#beta <- drawSamples(n, h = 5e-3)
beta <- drawSamples(n, h = 5e-3)

beta.normal <- t(rmvnorm(n/3, mean = rowMeans(beta), sigma = diag(1e3, 3, 3)))

beta <- cbind(beta, beta.normal)

beta.normal <- t(rmvnorm(n/3, mean = rowMeans(beta), sigma = diag(1, 3, 3)))

beta <- cbind(beta, beta.normal)

beta.normal <- t(rmvnorm(n/3, mean = rowMeans(beta), sigma = diag(1e-2, 3, 3)))

beta <- cbind(beta, beta.normal)

#beta <- t(rmvnorm(n, mean = numeric(3), sigma = diag(1, 3, 3)))
for (i in (n + 1):(2*n)) {
  print(i)
  log.F.samples[[toString(beta[,i])]] <- 0
  for (j in 1:length(Y)) {
    log.F.samples[[toString(beta[,i])]] <- log.F.samples[[toString(beta[,i])]] + as.numeric(Y[j]*t(Dx[j,]) %*% beta[,i] - exp(t(Dx[j,]) %*% beta[,i]))
  }
}

for (i in (n + 1):(2*n)) {
  print(i)
  H.samples[[toString(beta[,i])]] <- matrix(0, length(Dx[1,]), length(Dx[1,]))
  for (j in 1:length(Y)) {
    H.samples[[toString(beta[,i])]] <- H.samples[[toString(beta[,i])]] + as.numeric(exp(t(Dx[j,]) %*% beta[,i])) * Dx[j,] %*% t(Dx[j,])
  }
}

approx <- function(b, k = 5) {
  D <- list()
  # for (i in 1:length(beta[1,])) {
  #   D[[toString(as.vector(beta[,i]))]] <- norm(beta[,i] - b, type = "2")
  # }
  B <- names(log.F.samples)
  for (i in 1:length(B)) {
    be <- as.numeric(str_split(B[i], ", ")[[1]])
    D[[toString(as.vector(be))]] <- norm(be - b, type = "2")
  }
  sorted_keys <- names(sort(unlist(D)))[1:k]
  sorted_values <- sort(unlist(D))[1:k]
  s <- 0
  h <- matrix(0, length(Dx[1,]), length(Dx[1,]))
  d <- 0
  for (i in 1:k) {
    s <- s + log.F.samples[[sorted_keys[i]]]/(sorted_values[i])^2
    h <- h + H.samples[[sorted_keys[i]]]/(sorted_values[i])^2
    d <- d + 1/(sorted_values[i])^2
  }
  L <- list()
  L[["log.f"]] <- s/d
  L[["H"]] <- h/d
  return(L)
}


# sorted_keys <- names(sort(-unlist(log.F.samples)))[1:k]
# sorted_values <- -sort(-unlist(log.F.samples))[1:k]
# sorted_keys
# g$coefficients

be <- as.vector(rmvnorm(1, c(1, 4, 2), diag(3, 3, 3)))
be <- c(-0.5, 1.2, -0.8)
log.f.real(be)
as.numeric(approx(be, k = 150)$log.f)
#H.samples
#save(beta, log.F.samples, H.samples, file = "Model Approx Random.Rdata")
save(beta, log.F.samples, H.samples, file = "Model Approx to the Fourth Normal Random.Rdata")
