###############################################################################
# update_psi.R
#
# Updates the spike-and-slab allocation indicators for the chromosome
# arm-specific inclusion probabilities.
#
# For each chromosome arm g, psi_g[g] indicates whether the corresponding
# inclusion probability beta_g[g] is assigned to the slab (psi_g = 1) or
# the spike (psi_g = 0) component of the prior. Given the current values
# of beta_g and the mixing probability pi, each psi_g[g] is sampled from
# its Bernoulli full conditional distribution.
#
# Inputs:
#   beta_g : A numeric vector of chromosome arm-specific inclusion
#            probabilities.
#   b1     : First shape parameter of the slab Beta prior.
#   b2     : Second shape parameter of the slab Beta prior.
#   c1     : First shape parameter of the spike Beta prior.
#   c2     : Second shape parameter of the spike Beta prior.
#   pi_s   : Mixing probability of the spike-and-slab prior.
#
# Output:
#   A binary vector indicating the updated spike-and-slab allocation for
#   each chromosome arm.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# Function to update all psi_g given vector of beta_g values and parameters
update_psi <- function(beta_g, b1, b2, c1, c2, pi_s) {
  
  # Precompute log normalizing constants (computed once, not inside loop)
  log_B_b <- lbeta(b1, b2)   # log B(b1, b2)
  log_B_c <- lbeta(c1, c2)   # log B(c1, c2)
  
  # Initialize a vector to store updated psi_g values
  psi <- numeric(length(beta_g))
  
  # Iterate over each beta_g to update psi_g
  for (g in 1:length(beta_g)) {  # fixed: was beta, should be beta_g
    
    # Compute log probabilities for psi_g = 1 and psi_g = 0
    # including Beta normalizing constants which do NOT cancel here
    log_prob_1 <- (b1 - 1) * log(beta_g[g]) + (b2 - 1) * log(1 - beta_g[g]) 
    + log(pi_s)-log_B_b
    log_prob_0 <- (c1 - 1) * log(beta_g[g]) + (c2 - 1) * log(1 - beta_g[g]) 
    + log(1 - pi_s) - log_B_c
    
    # Normalized probability for psi_g = 1 using log-sum-exp for numerical stability
    log_sum   <- sumLog(c(log_prob_0, log_prob_1))
    prob      <- exp(log_prob_1 - log_sum)
    
    # Draw from Bernoulli distribution to determine psi_g
    psi[g] <- rbinom(1, 1, prob)
  }
  return(psi)
}