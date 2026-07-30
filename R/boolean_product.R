###############################################################################
# boolean_product.R
#
# Computes the Boolean matrix product of two binary matrices.
#
# Given a binary loading matrix W (K × R) and a binary signature matrix
# H (R × G), the Boolean product Z = W ∘ H is defined by
#
#   Z[k, g] = OR_{r=1}^R (W[k, r] AND H[r, g]),
#
# where multiplication is replaced by the logical AND operator and
# addition is replaced by the logical OR operator.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################

boolean_product <- function(W, H) {
  
  ## Check inputs
  stopifnot(is.matrix(W))
  stopifnot(is.matrix(H))
  stopifnot(ncol(W) == nrow(H))
  
  K <- nrow(W)
  G <- ncol(H)
  R <- ncol(W)
  
  X <- matrix(0L, nrow = K, ncol = G)
  
  for (r in seq_len(R)) {
    X <- X | (W[, r, drop = FALSE] %*% H[r, , drop = FALSE])
  }
  
  storage.mode(X) <- "integer"
  
  return(X)
}