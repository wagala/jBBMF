###############################################################################
# compute_log_postW.R
#
# Computes the log-posterior probability of a binary loading matrix.
#
# Given an observed binary data matrix X, a loading matrix W, a shared
# pattern matrix H, the patient-specific activation probabilities alpha,
# and the observation model parameters, this function computes the
# log-posterior probability of W. The log-posterior is obtained by summing
# the log-likelihood of the observed data and the log-prior probability of
# the loading matrix:
#
#   log P(W | X, H, alpha, p11, p10)
#     = log P(X | W, H, p11, p10)
#       + log P(W | alpha),
#
# up to an additive normalizing constant.
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
#   A numeric scalar containing the log-posterior probability of W,
#   up to an additive constant.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# compute_log_postW function
compute_log_postW <- function(X,W, H, alpha, p11, p10) {
  loglikeW <- computeLoglikelW(X,W, H, p11, p10)
  logpriorW <- compute_log_priorW(W, alpha)
  return(loglikeW + logpriorW)
}