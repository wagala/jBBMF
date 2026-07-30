###############################################################################
# update_H_Block_columns.R
#
# Updates the shared binary pattern matrix H using column-wise block Gibbs
# sampling.
#
# For each chromosome arm g, this function jointly updates all R latent-factor
# indicators in the column H[, g]. It enumerates the 2^R possible binary
# configurations for that chromosome arm and evaluates the log-posterior
# probability of H under each configuration using compute_log_posteriorH().
#
# The resulting log-posterior values are normalized using the numerically
# stable log-sum-exp calculation implemented in sumLog(). One complete column
# configuration is then sampled from the corresponding categorical
# distribution.
#
# Because H is shared across the diagnosis and relapse datasets, each candidate
# column configuration is evaluated using the joint likelihood of X1 and X2,
# together with the Bernoulli prior determined by beta[g].
#
# This block update is most practical when the number of latent factors R is
# relatively small, since the number of possible column configurations grows
# as 2^R.
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

#-------------------------------------------------------------------------------
update_H_Block_columns<- function(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta) {
  R <- nrow(H)  # Number of rows
  G <- ncol(H)  # Number of columns
  
  for (g in 1:G) {
    # Generate all possible configurations for the current column
    configurations <- expand.grid(rep(list(c(0, 1)), R))
    
    log_posteriors <- apply(configurations, 1, function(col_config) {
      H[, g] <- as.numeric(col_config)
      # Hypothetical function to compute log posterior
      #compute_log_posteriorH(W, H, X, beta, p11, p10) #original
      compute_log_posteriorH(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta)
    })
    
    
    # # Normalize the probabilities using the sumLog function for stability
    # log_sum_exp <- (log_posteriors)
    # probabilities <- exp(log_posteriors - log_sum_exp)
    
    # Convert log posteriors to probabilities using log-sum-exp for stability
    log_post = sumLog(log_posteriors)
    probs = exp(log_posteriors -log_post)
    # Use rmultinom to sample a configuration based on computed probabilities
    sampled_config_index <- which(rmultinom(1, 1, probs) == 1)
    H[, g] <- as.numeric(configurations[sampled_config_index, ])
  }
  
  return(H)
}