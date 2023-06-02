#' Explode dimension of a 2d matrix to 3d array
#'
#' @description
#' Turns a 2d matrix *A* of size *s* x *c* into a 3d array by replicating
#' rows of the matrix *A* either vertically or horizontally.
#'
#' This function takes a matrix *A*
#' ```
#' A = [[a_11, ..., a_1c],
#'      [... , ..., ... ]
#'      [a_s1, ..., a_sc]]
#' ```
#'
#' transforms each row `r` to form a separate 2d matrix *Ar*
#' based on the `byrow` parameters, and returns a 3d array
#' `[A1, ..., As]`.
#'
#' @param A
#' a `matrix`.
#'
#' @param byrow
#' states if the rows of the matrix should be replicated
#' _vertically_ (`byrow=FALSE`)
#' or _horizontally_ (`byrow=TRUE`).
#' Note the default behaviour is to replicate _horizontally_.
#'
#' _vertical_ replication means each *Ar* looks like
#' ```
#' A_r = [[a_r1, ..., a_r1],
#'        [... , ..., ... ],
#'        [a_rc, ..., a_rc]],
#' ```
#' while _horizontal_ replication means each *Ar* looks like
#' ```
#' A_r = [[a_r1, ..., a_rc],
#'        [... , ..., ... ],
#'        [a_r1, ..., a_rc]].
#' ```
#'
#' @return a 3d `array` with appropriately replicated rows.
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4), ncol=2)
#' # > A
#' #      [,1] [,2]
#' # [1,]    1    3
#' # [2,]    2    4
#' B <- explode_dimension(A)
#' # > B
#' # , , 1
#' #
#' #      [,1] [,2]
#' # [1,]    1    1
#' # [2,]    3    3
#' #
#' # , , 2
#' #
#' #      [,1] [,2]
#' # [1,]    2    2
#' # [2,]    4    4
explode_dimension <- function(A, byrow=FALSE) {
  stopifnot(is.matrix(A))
  apply(A, 1,
        function(x) matrix(x, ncol(A), ncol(A), byrow=byrow), simplify=FALSE) |>
    unlist() |>
    array(dim=c(ncol(A), ncol(A), nrow(A)))
}

update_cluster_centers <-
  function(
    U,
    X,
    alpha=NULL,
    F_=NULL #,
    # h_indices=  [TODO: not needed? can calc from F_]
  ) {
    if (is.null(alpha)) {
      V <-
        t(sweep(
          t(X) %*% U^2,
          2,
          colSums(U^2),
          "/"
        ))
    } else {    # need to calculate $\Phi$
      UF <- alpha * (U-F_)^2
      i_indices <- which(rowSums(F_) != 0)
      j_indices <- 1:nrow(F_)
      h_indices <- setdiff(j_indices, i_indices)
      UF[h_indices,] <- 0.
      Phi <- U^2 + UF

      V <- t(sweep(
        t(X) %*% Phi,
        2,
        colSums(Phi),
        "/"
      ))
    }

    return(V)
  }

#' Calculating distances between two matrices
#'
#' A base R approach to calculate distances between
#' rows of two matrices A and B that are of the same ncol.
#'
#' TODO it should be moved to tests to compare with other, externally
#' imported methods because this function turned out to be very slow.
#'
#' @param X
#' a matrix of dimension (N, p)
#'
#' @param V
#' a matrix of dimension (C, p)
#'
#' @returns a matrix of dimension (N, C)
#'
calculate_distances <-
  function(
    X,
    V
  ) {
    process_distance <- function(x, y) {
      output <- as.matrix(dist(rbind(x, y)))
      output[2:nrow(output), 1]
    }

    D_ <- apply(X, 1, function(x, y) process_distance(x, y), V, simplify=FALSE)
    D <- t(do.call(cbind, D_))

    return(D)
  }

#' Estimated U matrix with memberships
#'
#' The trick is to use 3d arrays for efficient calculations.
#' Let's take a D^2 matrix
#' ```
#' D62 = [[5, 145],
#'        [25, 85],
#'        [61, 41]]
#' ```
#'
#' `D_nominator` will be a 3d array of dimensions (2, 2, 3):
#' it contains 3 matrices 2x2.
#' Let's take a single 2x2 matrix from the array `D_nominator[,,1]`
#' ```
#'      [,1] [,2]
#' [1,]    5    5
#' [2,]  145  145
#' ```
#' It is effectively
#' ```
#'           [,1]      [,2]
#' [1,]    d^2_11    d^2_11
#' [2,]    d^2_12    d^2_12
#' ```
#' Let's now take `D_denominator[,,1]`
#' ```
#'      [,1] [,2]
#' [1,]    5  145
#' [2,]    5  145
#' ```
#' It is effectively
#' ```
#'           [,1]      [,2]
#' [1,]    d^2_11    d^2_12
#' [2,]    d^2_11    d^2_12
#' ```
#' If we now divide `D_nominator[,,1]` by `D_denominator[,,1]`, we get
#' ```
#'          [,1]     [,2]
#' [1,]    5/5    5/145
#' [2,]    145/5  145/145
#' ```
#' which is
#'                  [,1]             [,2]
#' [1,]    d^2_11/d^2_11    d^2_11/d^2_12
#' [2,]    d^2_12/d^2_11    d^2_12/d^2_12
#' The `rowSums` of the above will give us a vector
#' ```
#' [1] e_11 e_12
#' ```
#' *almost* the first row of evidence matrix, which is the inverse of the above
#' i.e. a vector `[1] 1/e_11 1/e_12`.
#'
#' The vectorized operations on entire `D_nominator` and `D_denominator`
#' follow the above logic.
#'
#'
update_memberships <-
  function(
    X,
    V,
    F_,
    alpha,
    fun.distances
  ) {
    D <- fun.distances(X, V)^2
    # TODO we need a thorough mathematical description here
    #
    D_nominator <- explode_dimension(D)
    D_denominator <- explode_dimension(D, byrow=TRUE)

    E_reciprocal <- t(apply(D_nominator/D_denominator, c(3), rowSums))
    E <- 1/E_reciprocal

    if (is.null(alpha)) {
      return(E)
    } else {
      i_indices <- which(rowSums(F_) != 0)
      M <- matrix(1, nrow(F_), ncol(F_))
      M[i_indices, ] <- 1/(1+alpha)

      ALB = F_*(alpha/(1+alpha))

      return(M*E + ALB)
    }
  }

#' Semi-Supervised Fuzzy C-Means model.
#'
#' @description
#' If *alpha* and *F_* are not supplied (their default values are `NULL`),
#' then a regular unsupervised Fuzzy C-Means algorithm is fitted.
#'
#' @param X
#' a matrix *X* with predictor variables.
#'
#' @param C
#' a number of clusters to find.
#'
#' @param U
#' optionally: a first memberships matrix to initialize the algorithm.
#' Used mainly for reproducibility to compare calculations with other packages
#' (e.g. in Python).
#'
#' @param fun.distances
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @export
#'
#' @examples
#' # simulate 50 obs from N2((5, 5)) and 50 obs from N2((7, 7))
#' library(MASS)
#' X <- rbind(
#'   MASS::mvrnorm(50, mu=c(5, 8), Sigma=matrix(c(3, 0, 0, 3), ncol=2)),
#'   MASS::mvrnorm(50, mu=c(7, 10), Sigma=matrix(c(3, 0, 0, 3), ncol=2))
#'   )
#'
#' # simulate supervision for 10% of each class
#' F_ <- matrix(0, nrow=100, ncol=2)
#' F_[sample(1:50, 10), 1] <- 1
#' F_[sample(51:100, 10), 2] <- 1
#'
#' model <- SSFCM(X=X, C=2)
#' model.ss <- SSFCM(X=X, C=2, alpha=1, F_=F_)
#'
#' acc.unsupervised.1 <- sum(apply(model$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#' acc.unsupervised.2 <- sum(apply(model$U, 1, which.max) == c(rep(2, 50), rep(1, 50)))
#' acc.supervised <- sum(apply(model.ss$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#'
SSFCM <- function(
    X,
    C,
    U=NULL,
    max_iter=200,
    conv_criterion=1e-4,
    fun.distances=rdist::cdist,
    alpha=NULL,
    F_=NULL
) {

  # random U if not supplied
  if (is.null(U)) {
    U <- matrix(runif(nrow(X)*C), ncol=C)
  }

  # normalize U
  U <- t(apply(U, 1, function(x) x / sum(x)))

  # counter
  counter = 0

  # calculations loop
  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U

    V <- update_cluster_centers(
      U=U_previous_iter,
      X=X,
      alpha=alpha,
      F_=F_)

    U <- update_memberships(
      X=X,
      V=V,
      F_=F_,
      alpha=alpha,
      fun.distances=fun.distances)

    conv_iter <- base::norm(U - U_previous_iter, type="F")

    if (conv_iter < conv_criterion) {
      break
    }
  }
  return(list(U=U, counter=counter))
}
