###############################################################################
# update_H.R
#
# Updates the shared binary pattern matrix H using element-wise Gibbs
# sampling.
#
# For each latent factor r and chromosome arm g, this function samples the
# indicator H[r, g] from its Bernoulli full conditional distribution. The
# conditional probability is obtained by evaluating the log-posterior of H
# under the two possible values:
#
#   H[r, g] = 0
#   H[r, g] = 1.
#
# Each candidate value is assessed using the joint diagnosis–relapse
# likelihood together with the Bernoulli prior determined by beta[g].
# The two log-posterior values are normalized using the numerically stable
# log-sum-exp calculation implemented in sumLog(). The updated value of
# H[r, g] is then sampled from the resulting Bernoulli distribution.
#
# Inputs:
#   X1     : A K × G binary diagnosis data matrix.
#   X2     : A K × G binary relapse data matrix.
#   W1     : A K × R binary diagnosis loading matrix.
#   W2     : A K × R binary relapse loading matrix.
#   H      : An R × G binary shared pattern matrix.
#   p1_11  : Probability of observing a 1 in X1 given X1_tilde = 1.
#   p1_10  : Probability of observing a 1 in X1 given X1_tilde = 0.
#   p2_11  : Probability of observing a 1 in X2 given X2_tilde = 1.
#   p2_10  : Probability of observing a 1 in X2 given X2_tilde = 0.
#   beta   : A numeric vector of length G containing the chromosome arm-
#            specific inclusion probabilities.
#
# Output:
#   The updated R × G binary shared pattern matrix H.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# update_H function
update_H <- function(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta) {
  G <- ncol(H)
  R <- nrow(H)
  for (r in 1:R) {
    for (g in 1:G) {
      original_h_rg <- H[r, g]
      
      # Compute the log posterior when h_rg = 0
      H[r, g] <- 0
      log_posterior_0 <- compute_log_posteriorH(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta)
      
      # Compute the log posterior when h_rg = 1
      H[r, g] <- 1
      log_posterior_1 <- compute_log_posteriorH(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta)
      
      # Reset H to its original value
      H[r, g] <- original_h_rg
      
      # Compute the normalized probabilities
      log_propH <- log_posterior_1 - sumLog(c(log_posterior_0, log_posterior_1))
      
      # Sample H
      p_H <- exp(log_propH)
      H[r, g] <- rbinom(1, 1, p_H)
    }
  }
  return(H)
}