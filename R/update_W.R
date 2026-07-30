###############################################################################
# update_W.R
#
# Updates the binary loading matrix W using element-wise Gibbs sampling.
#
# For each patient k and latent factor r, this function samples the loading
# indicator W[k, r] from its Bernoulli full conditional distribution. The
# conditional probability is obtained by evaluating the log-posterior of W
# under the two possible values:
#
#   W[k, r] = 0
#   W[k, r] = 1.
#
# The two log-posterior values are normalized using the numerically stable
# log-sum-exp calculation implemented in sumLog(). The updated value of
# W[k, r] is then sampled from the resulting Bernoulli distribution.
#
# This function is generic and can be used to update either the diagnosis
# loading matrix W1 or the relapse loading matrix W2 when the corresponding
# data matrix, activation probabilities, and observation parameters are
# supplied.
#
# Inputs:
#   X     : A K × G observed binary data matrix.
#   W     : A K × R binary loading matrix.
#   H     : An R × G binary shared pattern matrix.
#   alpha : A numeric vector of length K containing the patient-specific
#           activation probabilities.
#   p11   : Probability of observing a 1 given X_tilde = 1.
#   p10   : Probability of observing a 1 given X_tilde = 0.
#
# Output:
#   The updated K × R binary loading matrix W.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#update_W function
update_W <- function( X,W, H, alpha, p11, p10) {
  K <- nrow(W)
  R <- ncol(W)
  for (k in 1:K) {
    for (r in 1:R) {
      original_w_kr <- W[k, r]
      
      # Compute the log posterior when w_kr = 0
      W[k, r] <- 0
      log_posterior_0 <- compute_log_postW( X,W, H, alpha, p11, p10)
      
      # Compute the log posterior when w_kr = 1
      W[k, r] <- 1
      log_posterior_1 <- compute_log_postW(X, W, H,  alpha, p11, p10)
      
      # Reset W to its original value
      W[k, r] <- original_w_kr
      
      # Compute the normalized probabilities
      log_propW <- log_posterior_1 - sumLog(c(log_posterior_1, log_posterior_0))
      p_W <- exp(log_propW) # Normalized Probability
      
      # Update W based on the normalized probabilities
      W[k, r] <- rbinom(1, 1, p_W)
    }
  }
  return(W)
}
#update_W function
update_W <- function( X,W, H, alpha, p11, p10) {
  K <- nrow(W)
  R <- ncol(W)
  for (k in 1:K) {
    for (r in 1:R) {
      original_w_kr <- W[k, r]
      
      # Compute the log posterior when w_kr = 0
      W[k, r] <- 0
      log_posterior_0 <- compute_log_postW( X,W, H, alpha, p11, p10)
      
      # Compute the log posterior when w_kr = 1
      W[k, r] <- 1
      log_posterior_1 <- compute_log_postW(X, W, H,  alpha, p11, p10)
      
      # Reset W to its original value
      W[k, r] <- original_w_kr
      
      # Compute the normalized probabilities
      log_propW <- log_posterior_1 - sumLog(c(log_posterior_1, log_posterior_0))
      p_W <- exp(log_propW) # Normalized Probability
      
      # Update W based on the normalized probabilities
      W[k, r] <- rbinom(1, 1, p_W)
    }
  }
  return(W)
}