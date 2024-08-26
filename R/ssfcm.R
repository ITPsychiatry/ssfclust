#' Creates DHE (stands for "distances horizontally exploded") and DVE
#' (stands for "distances vertically exploded") matrices.
#' See vignette `vectorized-equations` for more details.
#'
#' @param A Matrix of size N x c.
#' @param vertical Boolean switch.
#' If `TRUE`, create DVE (vertical explosion).
#' If `FALSE`, create DHE (horizontal explosion).
#'
#' @return Matrix of size Nc x c
#' @export
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4, 5, 6), ncol=2, byrow=TRUE)
#' AVE <- dheve(A, vertical=TRUE)
#' AHE <- dheve(A, vertical=FALSE)
#'
dheve <- function(A, vertical) {
  if (vertical == TRUE) {
    elements <- A[rep(1:nrow(A), each=ncol(A)), ]
  } else {
    elements <- matrix(c(t(A)))[, rep(1, ncol(A))]
  }
  return(elements)
}


#' Aggregates elements of DHE and DVE matrices in a step to build
#' evidence matrix E.
#' See vignette `vectorized-equations` for details.
#'
#' @param dhe DHE matrix of size Nc x c.
#' @param dve DVE matrix of size Nc x c.
#'
#' @return Matrix of size Nc x 1.
#' @export
#'
gamma <- function(dhe, dve) {
  1 / ((dhe/dve) %*% matrix(rep(1, ncol(dhe))))
}


#' Rearranges elements of input matrix from a block matrix with vertical blocks
#' (column vectors) to a block matrix with horizontal blocks (row vectors).
#' See vignette `vectorized-equations` for details.
#'
#' @param A Matrix of size Nc x 1.
#' @param c Number of columns in the wanted matrix.
#' Associated with the number of clusters.
#'
#' @return Matrix of size N x c.
#' @export
#'
phi <- function(A, c) {
  matrix(A, ncol=c, byrow=TRUE)
}


#' Calculates data evidence matrix E from distances matrix D.
#'
#' @param D Distances matrix of size N x c.
#'
#' @return Matrix of size N x c.
#' @export
#'
calculate_evidence <- function(D) {
  dve <- dheve(D, vertical=TRUE)
  dhe <- dheve(D, vertical=FALSE)
  phi(A=gamma(dhe, dve), c=ncol(dhe))
}


#' Estimated U matrix with memberships.
#'
#' @param X
#' a matrix *X* of dimension (N, p) containing predictor variables.
#'
#' @param V
#' a prototypes matrix of dimension (c, p)
#'
#' @param F_
#' the supervision  binary matrix of the same dimension as *U*.
#'
#' @param alpha
#' the scaling factor, a floating point > 0.
#'
#' @param function_dist
#' A function of two arguments: matrices X and V of the same
#' number of columns.
#' It should return a matrix of (nrow(X) x nrow(V)) of distances
#' between each row of X and all rows of V.
#' In case of Euclidean distance, the result should not be squared!
#'
estimate_U <-
  function(
    X,
    V,
    F_,
    alpha,
    function_dist,
    i_indices
  ) {
    D <- function_dist(X, V)^2
    E <- calculate_evidence(D)

    if (is.null(alpha)) {
      return(E)
    } else {
      M <- matrix((1 / (1 + alpha)), ncol = 1)[, rep(1, ncol(F_))]
      F_alpha = F_*(alpha/(1+alpha))

      return(M*E + F_alpha)
    }
  }


estimate_U_ALT <-
  function(
    X,
    V,
    F_,
    alpha,
    function_dist,
    i_indices
  ) {
    A <- matrix(alpha, ncol = 1)[, rep(1, ncol(F_))]

    D <- function_dist(X, V)^2
    E <- calculate_evidence(D+A)

    if (is.null(alpha)) {
      return(E)
    } else {
      superF_dummy <- create_dummy_superF(F_)
      DSI <- matrix(get_supervised(D, superF_dummy),
                    ncol = 1)[, rep(1, ncol(F_))]

      F_alpha <- A / (DSI + A)
      M <- 1 - F_alpha

      return(M*E + F_alpha * F_)
    }
  }


get_supervised <- function(dataM, supervisionM) {
  mask <- which(supervisionM != 0, arr.ind=TRUE)
  mask.sorted <- mask[order(mask[, "row"], decreasing=FALSE), ]
  return(dataM[mask.sorted])
}


create_dummy_superF <- function(superF) {
  superF_dummy <- superF
  h_indices <- which(rowSums(superF_dummy) == 0)
  dummy_fill <- c(rep(0, ncol(superF) - 1), 1)
  superF_dummy[h_indices, ] <- dummy_fill

  return(superF_dummy)
}


#' Indices (row, column) of matrix cells other than `i,s(i)`.
#'
#' @param supervisionM F matrix
#'
#' @return dataframe with `row, col` columns.
#' @export
#'
#' @examples
#' F_ <- matrix(c(0, 0, 0, 1, 1, 0), nrow = 3, byrow = TRUE)
#' indices_to_zero_out <- get_non_isi_matrix_cells(B)
#'
get_non_isi_matrix_cells <- function(supervisionM) {
  mask <- which(supervisionM == 0, arr.ind = TRUE)
  mask_sorted <- mask[order(mask[, "row"], decreasing = FALSE), ]
  return(mask_sorted)
}




#' Equation to calculate clusters' prototypes matrix $\hat{V}$.
#'
#' @param Phi Matrix with weights of size N x c.
#'
#' @param X Matrix with predictors of size N x p.
#'
#' @return Clusters' prototypes matrix of size c x p.
#' @export
#'
estimate_V <- function(Phi, X) {
  Phi_tilde <- sweep(Phi, 2, colSums(Phi), "/")
  return(t(t(X) %*% Phi_tilde))
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
#' @param function_dist
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
#' model_ss <- SSFCM(X=X, C=2, alpha=1, F_=F_)
#'
#' acc_unsupervised_1 <- sum(apply(model$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#' acc_unsupervised_2 <- sum(apply(model$U, 1, which.max) == c(rep(2, 50), rep(1, 50)))
#' acc_supervised <- sum(apply(model_ss$U, 1, which.max) == c(rep(1, 50), rep(2, 50)))
#'
SSFCM <- function(
    X,
    C,
    U = NULL,
    max_iter = 200,
    conv_criterion = 1e-4,
    function_dist = rdist::cdist,
    alpha = NULL,
    F_ = NULL,
    prototypes = TRUE,
    alt = FALSE
) {
  if (is.null(U)) {
    U <- matrix(runif(nrow(X)*C), ncol=C)
  }

  # Rows of U should sum up to 1
  U <- t(apply(U, 1, function(x) x / sum(x)))

  # Calculate indices once instead in each loop
  if (is.null(alpha)) {
    i_indices <- NA
  } else {
    i_indices <- which(rowSums(F_) != 0)
    h_indices <- which(rowSums(F_) == 0)
  }

  counter = 0
  U_history <- list()
  V_history <- list()
  Phi_history <- list()

  for (iter in 1:max_iter) {
    counter <- counter + 1
    U_previous_iter <- U

    Phi <- U_previous_iter^2

    # Modify `Phi` if running semi-supervised FCM
    if (!is.null(alpha) && prototypes) {
      U_alpha <- alpha * (U_previous_iter - F_)^2
      U_alpha[h_indices, ] <- 0
      Phi <- Phi + U_alpha
    }

    Phi_history[[counter]] <- Phi

    V <- estimate_V(Phi, X)

    V_history[[counter]] <- V

    if (alt) {
      U <- estimate_U_ALT(
        X=X,
        V=V,
        F_=F_,
        alpha=alpha,
        function_dist=function_dist,
        i_indices=i_indices)
    } else {
      U <- estimate_U(
        X=X,
        V=V,
        F_=F_,
        alpha=alpha,
        function_dist=function_dist,
        i_indices=i_indices)
    }

    U_history[[counter]] <- U

    conv_iter <- base::norm(U - U_previous_iter, type="F")

    if (conv_iter < conv_criterion) {
      break
    }
  }

  z <- list(
    U = U,
    V = V,
    function_dist = function_dist,
    counter = counter,
    V_history = V_history,
    U_history = U_history,
    Phi_history = Phi_history
  )

  class(z) <- "ssfcm"

  return(z)
}


#' @export
estimate_T <-
  function(
    X ,
    V,
    F_,
    alpha,
    function_dist,
    i_indices,
    gammas  # vector
  ) {
    D <- function_dist(X, V)^2
    G <- matrix(gammas, nrow = 1)[rep(1, nrow(F_)), ]
    M1 <- F_ * alpha
    M2 <- matrix(rowSums(F_))[, rep(1, ncol(F_))] * alpha + 1
    Tm <- (G + (M1 * D) )/ (G + (M2 * D))

    return(Tm)
  }



#' @export
SSPCM <- function(
    X,
    C,
    U = NULL,
    gammas = NULL,
    max_iter = 200,
    conv_criterion = 1e-4,
    function_dist = rdist::cdist,
    alpha = NULL,
    F_ = NULL
) {
  if (is.null(U)) {
    Tm <- matrix(runif(nrow(X)*C), ncol=C)
  } else{
    Tm <- U
  }

  # Rows of U should sum up to 1
  Tm <- t(apply(Tm, 1, function(x) x / sum(x)))

  if (is.null(gammas)) {
    gammas <- rep(1, C)
  }

  # Calculate indices once instead in each loop
  if (is.null(alpha)) {
    i_indices <- NA
  } else {
    i_indices <- which(rowSums(F_) != 0)
    h_indices <- which(rowSums(F_) == 0)
  }

  counter = 0
  U_history <- list()
  V_history <- list()
  Phi_history <- list()

  for (iter in 1:max_iter) {
    counter <- counter + 1
    Tm_previous_iter <- Tm

    Phi <- Tm_previous_iter^2

    # Modify `Phi` if running semi-supervised PCM
    if (!is.null(alpha)) {
      Tm_alpha <- alpha * (Tm_previous_iter - F_)^2
      Tm_alpha[h_indices, ] <- 0
      Phi <- Phi + Tm_alpha
    }

    Phi_history[[counter]] <- Phi

    V <- estimate_V(Phi, X)

    V_history[[counter]] <- V

    Tm <- estimate_T(
      X = X,
      V = V,
      F_ = F_,
      alpha = alpha,
      function_dist = function_dist,
      i_indices = i_indices,
      gammas = gammas
    )

    U_history[[counter]] <- Tm

    conv_iter <- base::norm(Tm - Tm_previous_iter, type="F")

    if (conv_iter < conv_criterion) {
      break
    }
  }

  z <- list(
    Tm = Tm,
    V = V,
    function_dist = function_dist,
    counter = counter,
    gammas = gammas,
    V_history = V_history,
    U_history = U_history,
    Phi_history = Phi_history
  )

  class(z) <- "sspcm"

  return(z)
}


#' Title
#'
#' @param object
#' @param newdata
#'
#' @return
#'
#' @export
#' @examples
predict.ssfcm <- function(object, newdata) {
  output <- estimate_U(
    X = newdata,
    V = object$V,
    F_ = NULL,
    alpha = NULL,
    function_dist = object$function_dist,
    i_indices = NULL
  )
  return(output)
}


#' Title
#'
#' @param object
#' @param newdata
#'
#' @return
#'
#' @export
#' @examples
predict.sspcm <- function(object, newdata) {
  output <- estimate_Tm(
    X = newdata,
    V = object$V,
    F_ = NULL,
    alpha = NULL,
    function_dist = object$function_dist,
    i_indices = NULL,
    gammas = object$gammas
  )
  return(output)
}
