library(kernlab)
library(rdist)

my_kernel <- function(a, b, sigma) {
  euclid_norm <- sqrt(sum((a - b)^2))^2
  output <- exp(-sigma*euclid_norm)
  return(output)
}

X <- rnorm(12) |> matrix(ncol = 2)
V <- rnorm(6, mean = 1) |> matrix(ncol = 2)

my_kernel(X[2, ], V[1, ], sigma = 1)
kernlab::rbfdot(sigma = 1)(X[2, ], V[1, ])



# checking

function_dist=rdist::cdist

foo_outer <- function(x, function_dist) {
  output <- function_dist(x)
  return(output)
}

partial <- function(f, ...) {
  fixed_args <- list(...)
  function(...) {
    args <- c(fixed_args, list(...))
    do.call(f, args)
  }
}

x <- X[2, ]

foo_outer(cbind(X[2, ], V[1, ]), rdist::rdist)
foo_outer(cbind(X[2, ], V[1, ]), partial(rdist::rdist, metric = "manhattan"))

rdist::rdist(cbind(X[2, ], V[1, ]), metric = "manhattan")

rdist::rdist(cbind(X[2, ], V[1, ]), metric = partial(my_kernel, sigma = 1))
rdist::rdist(cbind(X[2, ], V[1, ]), metric = rbfdot(sigma = 1))






library(MASS)
X <- rbind(
  MASS::mvrnorm(50, mu=c(5, 8), Sigma=matrix(c(3, 0, 0, 3), ncol=2)),
  MASS::mvrnorm(50, mu=c(7, 10), Sigma=matrix(c(3, 0, 0, 3), ncol=2))
)

# simulate supervision for 10% of each class
F_ <- matrix(0, nrow=100, ncol=2)
F_[sample(1:50, 10), 1] <- 1
F_[sample(51:100, 10), 2] <- 1

model_ss <- SSFCM(X=X, C=2, alpha=1*rowSums(F_), F_=F_,
                  function_dist = partial(rdist::cdist, metric = rbfdot(sigma = 0.1)),
                  kernel = TRUE)

model_ss2 <- SSFCM(X=X, C=2, alpha=1*rowSums(F_), F_=F_,
                  function_dist = partial(rdist::cdist, metric = rbfdot(sigma = 0.05)),
                  kernel = TRUE)

model_ss3 <- SSFCM(X=X, C=2, alpha=1*rowSums(F_), F_=F_)

acc_supervised <- sum(apply(model_ss$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
acc_supervised2 <- sum(apply(model_ss2$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
acc_supervised3 <- sum(apply(model_ss3$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))


plot(model_ss$U, xlim = c(0, 1), ylim = c(0, 1))
points(model_ss2$U, xlim = c(0, 1), ylim = c(0, 1), col = "green")

data.frame(
  f = rowSums(F_),
  u = rowSums(model_ss2$U)
)




