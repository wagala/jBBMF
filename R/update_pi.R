###############################################################################
# update_pi.R
#
# Updates the mixing probability of the spike-and-slab prior for the shared
# pattern matrix H.
#
# The parameter pi represents the prior probability that a chromosome arm is
# assigned to the slab component of the spike-and-slab prior. Given the current
# allocation indicators psi_g, pi is sampled from its Beta full conditional
# distribution.
#
# Inputs:
#   psi_g : A binary vector indicating the prior component (spike or slab)
#           assigned to each chromosome arm.
#   d1    : First shape parameter of the Beta prior for pi.
#   d2    : Second shape parameter of the Beta prior for pi.
#
# Output:
#   A numeric scalar containing the updated mixing probability of the
#   spike-and-slab prior.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# Function to update \pi
update_pi <- function(psi_g, d1, d2) {
  # Calculate the total number of successes (psi_g = 1) and failures (psi_g = 0)
  sum_psi <- sum(psi_g)  # Total number where psi_g = 1
  G <- length(psi_g)     # Total number of psi_g values
  
  # Calculate the parameters for the Beta distribution
  alpha_post <- sum_psi + d1
  beta_post <- G - sum_psi + d2
  
  # Sample pi from the Beta distribution
  pi_s <- rbeta(1, alpha_post, beta_post)
  
  return(pi_s)
}