###############################################################################
# run_multiple_JBBMF.R
#
# Runs multiple Gibbs sampling chains for the joint Bayesian Boolean matrix
# factorization model.
#
# This function repeatedly calls runJBBMF() to generate independent MCMC
# chains for the paired diagnosis and relapse datasets. Each chain is assigned
# a distinct random-number seed to support reproducibility while preserving
# independence across chains.
#
# The chains are run in parallel using pbmclapply(), and an ETA-style progress
# indicator is displayed during computation. The returned object is a named
# list, with one element for each chain.
#
# Inputs:
#   X1       : A K × G binary diagnosis data matrix.
#   X2       : A K × G binary relapse data matrix.
#   W1_o     : An initial K × R binary diagnosis loading matrix.
#   W2_o     : An initial K × R binary relapse loading matrix.
#   H_o      : An initial R × G binary shared pattern matrix.
#   psi_o    : A binary vector of length G containing the initial spike-and-
#              slab allocation indicators.
#   n_iter   : Total number of Gibbs sampling iterations per chain.
#   burn_in  : Number of initial iterations discarded as burn-in.
#   a1, a2   : Shape parameters of the Beta prior for the patient-specific
#              activation probabilities.
#   b1, b2   : Shape parameters of the slab Beta distribution.
#   c1, c2   : Shape parameters of the spike Beta distribution.
#   d1, d2   : Shape parameters of the Beta prior for the spike-and-slab
#              mixing probability.
#   n_chains : Number of independent MCMC chains to run.
#
# Output:
#   A named list containing the output from runJBBMF() for each chain.
#   The list elements are labeled Chain_1, Chain_2, ..., Chain_n.
#
# Dependencies:
#   This function requires pbmclapply() and detectCores(), typically from the
#   pbmcapply and parallel packages, respectively.
#
# Notes:
#   Chain-specific seeds are set as 1001 + chain_id.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
##with a progress bar
run_multiple_JBBMF<- function(X1,X2, W1_o,W2_o, H_o, psi_o, n_iter, burn_in,
                              a1, a2, b1, b2, c1, c2, d1, d2,n_chains) {
  # Define a helper function to run a single chain with a specific seed
  run_chain <- function(chain_id) {
    set.seed(1001 + chain_id)  # Unique seed for each chain
    results <- runJBBMF(X1,X2, W1_o,W2_o, H_o, psi_o, n_iter, burn_in,
                        a1, a2, b1, b2, c1, c2, d1, d2)
    return(results)
  }
  
  # Use pbmclapply to run chains in parallel and show progress
  all_chain_results <- pbmclapply(1:n_chains, run_chain, mc.cores = detectCores(),
                                  mc.style = "ETA")
  
  # Rename list elements with chain labels
  names(all_chain_results) <- paste0("Chain_", 1:n_chains)
  
  return(all_chain_results)
}