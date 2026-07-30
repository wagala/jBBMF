###############################################################################
# computeLoglikelW.R
#
# Computes the log-likelihood of a single binary data matrix.
#
# Given an observed binary data matrix X, a loading matrix W, a shared
# pattern matrix H, and the observation model parameters, this function
# computes the log-likelihood of the observed data under the Boolean matrix
# factorization model.
#
# The latent Boolean signal X_tilde is first constructed as the Boolean
# product W ∘ H. The log-likelihood is then evaluated using the observation
# model
#
#   P(X = 1 | X_tilde = 1) = p11,
#   P(X = 1 | X_tilde = 0) = p10.
#
# This function evaluates the contribution of a single dataset (e.g.,
# diagnosis or relapse) to the overall model likelihood.
#
# Inputs:
#   X    : A K × G observed binary data matrix.
#   W    : A K × R binary loading matrix.
#   H    : An R × G binary shared pattern matrix.
#   p11  : Probability of observing a 1 given X_tilde = 1.
#   p10  : Probability of observing a 1 given X_tilde = 0.
#
# Output:
#   A numeric scalar containing the log-likelihood of the observed data.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#Compute the Total Log Likelihood
computeLoglikelW<- function(X, W, H,p11, p10){
  # Computing X tildes
  X_tilde <- as.matrix(W%&%(H))*1
  # Helper: log-likelihood 
  log_likhood  <- function(X, X_tilde, p11, p10) {
    loglike <-sum(
      X_tilde * X * log(p11) + X_tilde * (1 - X) * log(1 - p11) +
        (1 - X_tilde) * X * log(p10) + (1 - X_tilde) * (1 - X) * log(1 - p10)
    )
    return(loglike)
  }
  # Compute total log-likelihood
  LL <-log_likhood(X, X_tilde, p11, p10)
  return(LL)
}
