###############################################################################
# logpost_W2_W1.R
#
# Computes the log-posterior probability of the relapse loading matrix W2
# conditional on the diagnosis loading matrix W1.
#
# Given the relapse data matrix X2, the diagnosis loading matrix W1, the
# relapse loading matrix W2, the shared pattern matrix H, the observation
# model parameters, and the transition probabilities gamma11 and gamma00,
# this function computes the log-posterior probability of W2. The
# log-posterior is obtained by summing the log-likelihood of the relapse
# data and the conditional log-prior relating W2 to W1:
#
#   log P(W2 | X2, W1, H, gamma11, gamma00, p11, p10)
#     = log P(X2 | W2, H, p11, p10)
#       + log P(W2 | W1, gamma11, gamma00),
#
# up to an additive normalizing constant.
#
# Inputs:
#   W2      : A K × R binary relapse loading matrix.
#   W1      : A K × R binary diagnosis loading matrix.
#   H       : An R × G binary shared pattern matrix.
#   X2      : A K × G binary relapse data matrix.
#   p11     : Probability of observing a 1 given X2_tilde = 1.
#   p10     : Probability of observing a 1 given X2_tilde = 0.
#   gamma11 : Probability that W2 = 1 given W1 = 1.
#   gamma00 : Probability that W2 = 0 given W1 = 0.
#   seed    : Optional random seed (unused by this function).
#
# Output:
#   A numeric scalar containing the log-posterior probability of W2,
#   up to an additive normalizing constant.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#-------------------------------------------------------------------------------
logpost_W2_W1 <- function(W2, W1, H, X2, p11, p10,gamma11,gamma00,seed = NULL){
  ll <- computeLoglikelW(X2, W2, H,p11, p10)
  cpl <-log_W2_given_W1(W1, W2,gamma11,gamma00,seed = NULL, eps = 1e-12)
  log_post <- ll + cpl
  return(log_post)
  
}