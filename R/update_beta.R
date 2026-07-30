###############################################################################
# update_beta.R
#
# Updates the chromosome arm-specific inclusion probabilities for the shared
# pattern matrix H.
#
# For each chromosome arm g, beta[g] represents the probability that the arm
# belongs to a latent factor in the shared pattern matrix H. Given the current
# values of H, each beta[g] is sampled from its Beta full conditional
# distribution using the number of factors in which chromosome arm g is active
# or inactive.
#
# Inputs:
#   H  : An R × G binary shared pattern matrix.
#   b1 : First shape parameter of the Beta prior for beta.
#   b2 : Second shape parameter of the Beta prior for beta.
#
# Output:
#   A numeric vector of length G containing updated chromosome arm-specific
#   inclusion probabilities for H.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################

#-------------------------------------------------------------------------------
update_beta<-function(H,b1,b2){
  G<-dim(H)[2]
  beta <- rep(NA, G)
  #update beta
  for(g in 1:G){
    beta[g] <- rbeta(1,shape1=sum(H[,g])+b1, shape2=sum(1-H[,g])+b2)
  }
  return(beta)
}