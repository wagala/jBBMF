###############################################################################
# update_W_Block_rows.R
#
# Updates the binary loading matrix W using row-wise block Gibbs sampling.
#
# For each patient k, this function jointly updates all R loading indicators
# in the row W[k, ]. It first enumerates the 2^R possible binary factor
# configurations for that patient. For each configuration, the function
# evaluates the corresponding log-posterior probability using
# compute_log_postW().
#
# The resulting log-posterior values are normalized using the numerically
# stable log-sum-exp calculation implemented in sumLog(). One complete row
# configuration is then sampled from the resulting categorical distribution.
#
# This block update allows the factor indicators for a patient to be sampled
# jointly rather than one element at a time. It is therefore most practical
# when the number of latent factors R is relatively small, since the number
# of possible row configurations increases as 2^R.
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
update_W_Block_rows <- function(X, W, H, alpha, p11, p10) {
  K <- nrow(W)  # Number of rows
  R <- ncol(W)  # Number of columns
  
  for (k in 1:K) {
    # Generate all possible configurations for the current row
    configurations <- expand.grid(rep(list(c(0, 1)), R))
    
    log_posteriors <- apply(configurations, 1, function(row_config) {
      W[k, ] <- as.numeric(row_config)
      compute_log_postW(X, W, H, alpha, p11, p10)
    })
    
    # Normalize the log posteriors to probabilities
    log_post = sumLog(log_posteriors)
    probs = exp(log_posteriors - log_post)
    
    # Use a multinomial approach to sample a configuration based on computed probabilities
    sampled_config_index <- which(rmultinom(1, 1, probs) == 1)
    W[k, ] <- as.numeric(configurations[sampled_config_index, ])
  }
  
  return(W)
}