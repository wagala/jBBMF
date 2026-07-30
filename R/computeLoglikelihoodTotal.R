###############################################################################
# computeLoglikelihoodTotal.R
#
# Computes the joint log-likelihood of the diagnosis and relapse data.
#
# Given the observed data matrices X1 and X2, the corresponding loading
# matrices W1 and W2, the shared pattern matrix H, and the observation
# model parameters, this function computes the total log-likelihood of the
# JBBMF model.
#
# The latent Boolean signals X1_tilde and X2_tilde are first constructed as
# the Boolean products W1 ∘ H and W2 ∘ H, respectively. The log-likelihood
# contributions from the diagnosis and relapse datasets are then computed
# independently and summed to obtain the joint log-likelihood.
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
#
# Output:
#   A numeric scalar containing the joint log-likelihood of the diagnosis
#   and relapse data under the current model parameters.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#Compute the Total Log Likelihood
computeLoglikelihoodTotal<- function(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10){
  # Computing X tildes
  X1_tilde <- as.matrix(W1%&%(H))*1
  X2_tilde <- as.matrix(W2%&%(H))*1
  # Helper: log-likelihood 
  log_likhood  <- function(X, X_tilde, p11, p10) {
    loglike <-sum(
      X_tilde * X * log(p11) + X_tilde * (1 - X) * log(1 - p11) +
        (1 - X_tilde) * X * log(p10) + (1 - X_tilde) * (1 - X) * log(1 - p10)
    )
    return(loglike)
  }
  # Compute total log-likelihood
  LL1 <-log_likhood(X1, X1_tilde, p1_11, p1_10)
  LL2 <-log_likhood(X2, X2_tilde, p2_11, p2_10)
  
  total_log_likelihood <- sum(LL1 + LL2)
  return(total_log_likelihood)
}