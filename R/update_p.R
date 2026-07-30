###############################################################################
# update_p.R
#
# Updates the observation model parameters for a binary data matrix.
#
# Given the current latent loading matrix W, shared pattern matrix H, and
# observed binary data matrix X, this function updates the observation
# probabilities p11 and p10. The latent Boolean signal X_tilde is first
# computed as the Boolean product of W and H. The observation probabilities
# are then sampled from their Beta full conditional distributions:
#
#   p11 = P(X = 1 | X_tilde = 1),
#   p10 = P(X = 1 | X_tilde = 0).
#
# Inputs:
#   W : A K × R binary loading matrix.
#   H : An R × G binary shared pattern matrix.
#   X : A K × G observed binary data matrix.
#
# Output:
#   A numeric vector containing the updated observation probabilities
#   c(p11, p10).
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#-------------------------------------------------------------------------------
update_p <- function(W, H, X) {
  K <- nrow(W)
  G <- ncol(H)
  
  # Compute x_tilde using matrix multiplication and convert to binary
  X_tilde <- as.matrix(W %&% H)*1
  
  # Compute the values for p11 and p10 using vectorized operations
  m1 <- sum((1 - X_tilde) * X)
  m <- sum(1 - X_tilde)
  n1 <- sum(X_tilde * X)
  n <- sum(X_tilde)
  
  # Sample p11 and p10 from the Beta distribution
  p11 <- rbeta(1, shape1 = n1 + 1, shape2 = n - n1 + 1)
  p10 <- rbeta(1, shape1 = m1 + 1, shape2 = m - m1 + 1)
  return(c(p11, p10))
}