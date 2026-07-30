###############################################################################
# update_beta_g.R
#
# Updates the chromosome arm-specific inclusion probabilities for the shared
# pattern matrix H under the spike-and-slab prior.
#
# For each chromosome arm g, beta_g[g] represents the probability that the
# arm is active in the shared pattern matrix H. The corresponding indicator
# psi_g[g] determines which Beta prior is assigned to beta_g[g]. If psi_g[g]=1,
# beta_g[g] is sampled from the slab component; otherwise, it is sampled from
# the spike component. The posterior full conditional is obtained by combining
# the selected Beta prior with the current values of column g of H.
#
# Inputs:
#   H      : An R × G binary shared pattern matrix.
#   psi_g  : A binary vector indicating the prior component (spike or slab)
#            for each chromosome arm.
#   b1     : First shape parameter of the slab Beta prior.
#   b2     : Second shape parameter of the slab Beta prior.
#   c1     : First shape parameter of the spike Beta prior.
#   c2     : Second shape parameter of the spike Beta prior.
#
# Output:
#   A numeric vector of length G containing updated chromosome arm-specific
#   inclusion probabilities.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# Function to update beta_g vector
update_beta_g <- function(H, psi_g, b1, b2, c1, c2) {
  R <- nrow(H)  # Number of observations per group (number of rows)
  G <- ncol(H)  # Number of groups (number of columns)
  
  # Initialize a vector to store sampled beta_g values
  beta_g_samples <- numeric(G)
  # Calculate the number of successes for each group g
  H_g <- colSums(H)
  # Iterate over each group g
  for (g in 1:G) {
    if (psi_g[g] == 1) {
      # Sample from Beta(H_g + b1, R - H_g + b2) when psi_g = 1
      beta_g_samples[g] <- stats::rbeta(1, H_g[g] + b1, R - H_g[g] + b2)
    } else if (psi_g[g] == 0) {
      # Sample from Beta(H_g + c1, R - H_g + c2) when psi_g = 0
      beta_g_samples[g] <- stats::rbeta(1, H_g[g] + c1, R - H_g[g] + c2)
    } else {
      stop("psi_g must be either 0 or 1.")
    }
  }
  return(beta_g_samples)
}