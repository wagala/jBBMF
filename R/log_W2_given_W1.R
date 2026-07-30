###############################################################################
# log_W2_given_W1.R
#
# Computes the log-prior probability of the relapse loading matrix W2
# conditional on the diagnosis loading matrix W1.
#
# Given the diagnosis loading matrix W1, the relapse loading matrix W2,
# and the transition probabilities gamma11 and gamma00, this function
# computes
#
#   log P(W2 | W1, gamma11, gamma00),
#
# where
#
#   gamma11 = P(W2 = 1 | W1 = 1),
#   gamma00 = P(W2 = 0 | W1 = 0).
#
# The conditional log-prior is evaluated from the joint cell counts
# (n11, n10, n01, n00), where
#
#   n11 : W1 = 1 and W2 = 1,
#   n10 : W1 = 1 and W2 = 0,
#   n01 : W1 = 0 and W2 = 1,
#   n00 : W1 = 0 and W2 = 0.
#
# Inputs:
#   W1      : A K × R binary diagnosis loading matrix.
#   W2      : A K × R binary relapse loading matrix.
#   gamma11 : Probability that W2 = 1 given W1 = 1.
#   gamma00 : Probability that W2 = 0 given W1 = 0.
#   seed    : Optional random seed (unused by this function).
#   eps     : Small positive constant reserved for numerical stability
#             (unused by this function).
#
# Output:
#   A numeric scalar containing the conditional log-prior probability
#   log P(W2 | W1, gamma11, gamma00).
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# ------------------------------log P(W2 | W1, gamma11, gamma00) ---------------

log_W2_given_W1<- function(W1, W2,gamma11,gamma00, seed = NULL, eps = 1e-12) {
  cnt <- counts_for_W1_W2(W1, W2)
  g11 <- gamma11; g00 <-gamma00
  log_prob <- cnt["n11"] * log(g11) +
    cnt["n10"] * log1p(-g11) +
    cnt["n00"] * log(g00) +
    cnt["n01"] * log1p(-g00)
  return(unname(log_prob))
}