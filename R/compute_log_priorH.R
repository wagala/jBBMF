###############################################################################
# compute_log_priorH.R
#
# Computes the log-prior probability of the shared pattern matrix H.
#
# Given the binary shared pattern matrix H and the chromosome arm-specific
# inclusion probabilities beta, this function computes the log-prior
# probability of H under the independent Bernoulli prior
#
#   H_rg ~ Bernoulli(beta_g),
#
# where beta_g is the probability that chromosome arm g belongs to latent
# factor r.
#
# Inputs:
#   H    : An R × G binary shared pattern matrix.
#   beta : A numeric vector of length G containing the chromosome arm-
#          specific inclusion probabilities.
#
# Output:
#   A numeric scalar containing the log-prior probability of H.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#function to calculate the log prior for H
compute_log_priorH <- function(H, beta) {
  G <- ncol(H)
  log_beta <- log(beta)
  log_one_minus_beta <- log(1 - beta)
  logpriorH <- sum(H %*% log_beta + (1 - H) %*% log_one_minus_beta)
  return(logpriorH)
}