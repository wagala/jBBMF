###############################################################################
# gammas_W1_W2.R
#
# Updates the dependence parameters relating the diagnosis and relapse
# loading matrices.
#
# Given the diagnosis loading matrix W1 and the relapse loading matrix W2,
# this function samples the transition probabilities gamma11 and gamma00
# from their conjugate Beta posterior distributions. The updates are based
# on the joint cell counts (n11, n10, n01, n00) computed from W1 and W2:
#
#   gamma11 = P(W2 = 1 | W1 = 1),
#   gamma00 = P(W2 = 0 | W1 = 0).
#
# Assuming independent Beta hyperpriors,
#
#   gamma11 ~ Beta(u11, v11),
#   gamma00 ~ Beta(u00, v00),
#
# the posterior distributions are
#
#   gamma11 | W1, W2 ~ Beta(u11 + n11, v11 + n10),
#   gamma00 | W1, W2 ~ Beta(u00 + n00, v00 + n01).
#
# Inputs:
#   W1    : A K × R binary diagnosis loading matrix.
#   W2    : A K × R binary relapse loading matrix.
#   u11   : First shape parameter of the Beta prior for gamma11.
#   v11   : Second shape parameter of the Beta prior for gamma11.
#   u00   : First shape parameter of the Beta prior for gamma00.
#   v00   : Second shape parameter of the Beta prior for gamma00.
#   seed  : Optional random seed for reproducible sampling.
#   eps   : Small positive constant used to bound sampled probabilities
#           away from 0 and 1 for numerical stability.
#
# Output:
#   A named numeric vector containing:
#
#     gamma11 : Updated probability that W2 = 1 given W1 = 1.
#     gamma00 : Updated probability that W2 = 0 given W1 = 0.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# -------- Sample gammas from the *posterior* (conjugate Beta) -----------------
# Hyperpriors:  gamma11 ~ Beta(u11, v11),  gamma00 ~ Beta(u00, v00)
# Posterior:    gamma11 | W ~ Beta(u11+n11, v11+n10)
#               gamma00 | W ~ Beta(u00+n00, v00+n01)
gammas_W1_W2 <- function(W1, W2, u11, v11, u00, v00, seed = NULL, eps = 1e-12) {
  if (!is.null(seed)) set.seed(seed)
  cnt <- counts_for_W1_W2(W1, W2)
  a11 <- u11 + cnt["n11"]; b11 <- v11 + cnt["n10"]
  a00 <- u00 + cnt["n00"]; b00 <- v00 + cnt["n01"]
  g11 <- rbeta(1L, a11, b11)
  g00 <- rbeta(1L, a00, b00)
  g11 <- min(max(g11, eps), 1 - eps)   # clamp for log safety
  g00 <- min(max(g00, eps), 1 - eps)
  c(gamma11 = g11, gamma00 = g00)
}