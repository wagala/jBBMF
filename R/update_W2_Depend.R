###############################################################################
# update_W2_Depend.R
#
# Updates the relapse loading matrix W2 conditional on the diagnosis loading
# matrix W1 using element-wise Gibbs sampling.
#
# For each patient k and latent factor r, this function samples W2[k, r] from
# its Bernoulli full conditional distribution. The conditional probability is
# obtained by evaluating the log-posterior of W2 under the two possible values:
#
#   W2[k, r] = 0
#   W2[k, r] = 1.
#
# Each candidate value is assessed using the relapse-data likelihood together
# with the conditional prior linking W2 to W1 through gamma11 and gamma00:
#
#   gamma11 = P(W2 = 1 | W1 = 1),
#   gamma00 = P(W2 = 0 | W1 = 0).
#
# The two log-posterior values are normalized using the numerically stable
# log-sum-exp calculation implemented in sumLog(). The updated value of
# W2[k, r] is then sampled from the resulting Bernoulli distribution.
#
# Inputs:
#   W2      : A K × R binary relapse loading matrix.
#   W1      : A K × R binary diagnosis loading matrix.
#   H       : An R × G binary shared pattern matrix.
#   X2      : A K × G binary relapse data matrix.
#   p11     : Probability of observing a 1 in X2 given X2_tilde = 1.
#   p10     : Probability of observing a 1 in X2 given X2_tilde = 0.
#   gamma11 : Probability that W2 = 1 given W1 = 1.
#   gamma00 : Probability that W2 = 0 given W1 = 0.
#   seed    : Optional random seed passed to the conditional log-posterior
#             function.
#
# Output:
#   The updated K × R binary relapse loading matrix W2.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#update_W2 when X2 is dependent X1 
update_W2_Depend<- function(W2, W1, H, X2, p11, p10,gamma11,gamma00, seed = NULL) {
  K <- nrow(W2)
  R <- ncol(W2)
  for (k in 1:K) {
    for (r in 1:R) {
      original_w_kr <- W2[k, r]
      
      # Compute the log posterior when w_kr = 0
      W2[k, r] <- 0
      log_posterior_0 <- logpost_W2_W1(W2, W1, H, X2, p11, p10,gamma11,gamma00,seed = NULL)
      
      # Compute the log posterior when w_kr = 1
      W2[k, r] <- 1
      log_posterior_1 <- logpost_W2_W1(W2, W1, H, X2, p11, p10,gamma11,gamma00,seed = NULL)
      
      # Reset W to its original value
      W2[k, r] <- original_w_kr
      
      # Compute the normalized probabilities
      log_propW <- log_posterior_1 - sumLog(c(log_posterior_1, log_posterior_0))
      p_W <- exp(log_propW) # Normalized Probability
      
      # Update W based on the normalized probabilities
      W2[k, r] <- rbinom(1, 1, p_W)
    }
  }
  return(W2)
}