###############################################################################
# counts_for_W1_W2.R
#
# Computes the joint cell counts for the diagnosis and relapse loading
# matrices.
#
# Given the binary diagnosis loading matrix W1 and the binary relapse loading
# matrix W2, this function computes the four joint counts describing their
# element-wise agreement:
#
#   n11 : Number of entries where W1 = 1 and W2 = 1.
#   n10 : Number of entries where W1 = 1 and W2 = 0.
#   n01 : Number of entries where W1 = 0 and W2 = 1.
#   n00 : Number of entries where W1 = 0 and W2 = 0.
#
# These counts are sufficient statistics for updating the transition
# probabilities governing the dependence of the relapse loading matrix W2
# on the diagnosis loading matrix W1.
#
# Inputs:
#   W1 : A K × R binary diagnosis loading matrix.
#   W2 : A K × R binary relapse loading matrix.
#
# Output:
#   A named numeric vector containing the counts:
#
#     n11, n10, n01, n00.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# ------------------ Counts: n11, n10, n01, n00 --------------------------------
counts_for_W1_W2 <- function(W1, W2) {
  stopifnot(is.matrix(W1), is.matrix(W2), all(dim(W1) == dim(W2)))
  W1b <- (W1 != 0); W2b <- (W2 != 0)
  c(n11 = sum(W1b &  W2b),
    n10 = sum(W1b & !W2b),
    n01 = sum(!W1b &  W2b),
    n00 = sum(!W1b & !W2b))
}