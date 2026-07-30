###############################################################################
# compute_log_posteriorH.R
#
# Computes the log-posterior probability of the shared pattern matrix H.
#
# Given the diagnosis and relapse data matrices (X1 and X2), the corresponding
# loading matrices (W1 and W2), the shared pattern matrix H, the observation
# model parameters, and the chromosome arm-specific inclusion probabilities,
# this function computes the log-posterior probability of H. The log-posterior
# is obtained by combining the joint log-likelihood of the two datasets with
# the prior distribution on H:
#
#   log P(H | X1, X2, W1, W2, beta, p)
#     = log P(X1, X2 | W1, W2, H, p)
#       + log P(H | beta),
#
# up to an additive normalizing constant.
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
#   beta   : Chromosome arm-specific inclusion probabilities for H.
#
# Output:
#   A numeric scalar containing the log-posterior probability of H,
#   up to an additive normalizing constant.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#-------------------------------------------------------------------------------
#compute_log_posteriorH function
compute_log_posteriorH <- function(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta) {
  loglikeH<-computeLoglikelihoodTotal(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10)
  logpriorH <- compute_log_priorH(H, beta)
  return(loglikeH + logpriorH)
}