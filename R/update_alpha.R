###############################################################################
# update_alpha.R
#
# Updates the patient-specific activation probabilities for the diagnosis
# loading matrix W_1.
#
# For each patient k, alpha[k] represents the probability that a latent factor
# is active in the diagnosis loading matrix W_1. Given the current values of
# W_1, each alpha[k] is sampled from its Beta full conditional distribution
# using the number of active and inactive latent factors for patient k.
#
# Inputs:
#   W   : A K × R binary diagnosis loading matrix (W_1).
#   a1  : First shape parameter of the Beta prior for alpha.
#   a2  : Second shape parameter of the Beta prior for alpha.
#
# Output:
#   A numeric vector of length K containing updated patient-specific
#   activation probabilities for W_1.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#By row
update_alpha<-function(W,a1,a2){
  K<-dim(W)[1]
  alpha <- rep(NA, K)
  #update alpha
  for (k in 1:K) {
    alpha[k] <- rbeta(1, shape1=sum(W[k,])+a1, shape2= sum(1-W[k,])+a1)
  }
  return(alpha)
}