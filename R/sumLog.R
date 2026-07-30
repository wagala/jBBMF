###############################################################################
# sumLog.R
#
# Computes the logarithm of the sum of exponentials of a numeric vector using
# a numerically stable implementation of the log-sum-exp trick.
#
# This function is used to avoid numerical overflow and underflow when
# evaluating expressions of the form log(sum(exp(vec))), which frequently
# arise in Bayesian inference and probabilistic computations.
#
# Inputs:
#   vec : A numeric vector containing values on the log scale.
#
# Output:
#   A numeric scalar equal to log(sum(exp(vec))).
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################

# sumLog function
sumLog <- function(vec) {
  ord <- sort(vec, decreasing = TRUE)
  s <- ord[1]
  s <- s + sum(log1p(exp(ord[-1] - s)))
  return(s)
}