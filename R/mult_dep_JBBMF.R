###############################################################################
# mult_dep_JBBMF.R
#
# Runs multiple independent MCMC chains for the dependent Joint Bayesian
# Boolean Matrix Factorization model.
#
# This function fits the dependent JBBMF model across several parallel chains
# by repeatedly calling dependent_JBBMF(). Each chain uses the same observed
# data, initial matrices, prior hyperparameters, number of iterations, and
# burn-in period, but is assigned a distinct reproducible random seed.
#
# The chain-specific seed is defined as
#
#   chain_seed = seed + chain_id.
#
# Parallel computation is performed using pbmclapply(), with the number of
# available processor cores determined by detectCores(). The returned list is
# labeled by chain to facilitate convergence assessment and posterior
# combination.
#
# Inputs:
#   X1       : A K × G binary diagnosis data matrix.
#   X2       : A K × G binary relapse data matrix.
#   W1_o     : Initial K × R binary diagnosis loading matrix.
#   W2_o     : Initial K × R binary relapse loading matrix.
#   H_o      : Initial R × G binary shared pattern matrix.
#   psi_o    : Initial length-G vector of spike-and-slab indicators.
#   n_iter   : Total number of Gibbs sampling iterations per chain.
#   burn_in  : Number of initial iterations discarded as burn-in.
#   a1, a2   : Beta prior hyperparameters for the patient-specific activation
#              probabilities associated with W1.
#   b1, b2   : Beta prior hyperparameters for the slab component of the
#              chromosome arm-specific inclusion probabilities.
#   c1, c2   : Beta prior hyperparameters for the spike component of the
#              chromosome arm-specific inclusion probabilities.
#   d1, d2   : Beta prior hyperparameters for the spike-and-slab mixing
#              probability.
#   u11, v11 : Beta prior hyperparameters for gamma11.
#   u00, v00 : Beta prior hyperparameters for gamma00.
#   n_chains : Number of independent MCMC chains.
#   seed     : Base random seed used to construct chain-specific seeds.
#
# Output:
#   A named list containing the output from dependent_JBBMF() for each chain.
#   The list elements are labeled:
#
#     Chain_1, Chain_2, ..., Chain_n.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
mult_dep_JBBMF <- function(X1, X2,W1_o, W2_o, H_o, psi_o,n_iter,
                           burn_in,a1, a2, b1, b2, c1, c2, d1, d2,u11, 
                           v11, u00, v00,n_chains,seed) {
  
  # helper: one chain with unique, reproducible seed
  run_chain <- function(chain_id) {
    chain_seed <- seed + chain_id
    set.seed(chain_seed)# Unique seed for each chain
    results <- dependent_JBBMF(X1, X2,W1_o, W2_o, H_o, psi_o,n_iter, burn_in,a1, 
                               a2, b1, b2, c1, c2, d1, d2,u11, v11, u00, v00,seed = chain_seed)
    return(results)
  }
  #Use pbmclapply to run chains in parallel and show progress
  all_chain_results <- pbmclapply(1:n_chains, run_chain, mc.cores = detectCores(),
                                  mc.style = "ETA")
  # Rename list elements with chain labels
  names(all_chain_results) <- paste0("Chain_", 1:n_chains)
  
  return(all_chain_results)
}
#-------------------------------------------------------------------------------