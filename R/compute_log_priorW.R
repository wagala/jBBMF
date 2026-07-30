###############################################################################
# compute_log_priorW.R
#
# Computes the log-prior probability of a binary loading matrix.
#
# Given a binary loading matrix W and the corresponding vector of
# patient-specific activation probabilities alpha, this function computes the
# log-prior probability of W under the independent Bernoulli prior
#
#   W_kr ~ Bernoulli(alpha_k),
#
# where alpha_k is the probability that latent factor r is active for
# patient k.
#
# Inputs:
#   W     : A K × R binary loading matrix.
#   alpha : A numeric vector of length K containing the patient-specific
#           activation probabilities.
#
# Output:
#   A numeric scalar containing the log-prior probability of W.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# compute_log_priorW function
compute_log_priorW <- function(W, alpha) {
  logpriorW <- sum(W * log(alpha) + (1 - W) * log(1 - alpha))
  return(logpriorW)
}