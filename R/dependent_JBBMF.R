###############################################################################
# dependent_JBBMF.R
#
# Runs the Gibbs sampler for the dependent Joint Bayesian Boolean Matrix
# Factorization model.
#
# This function fits a two-view Boolean matrix factorization model in which
# the diagnosis and relapse data matrices, X1 and X2, are represented by
#
#   X1_tilde = W1 o H,
#   X2_tilde = W2 o H,
#
# where W1 and W2 are view-specific binary loading matrices and H is a shared
# binary pattern matrix. Dependence between the two loading matrices is
# introduced through
#
#   gamma11 = P(W2 = 1 | W1 = 1),
#   gamma00 = P(W2 = 0 | W1 = 0).
#
# At each MCMC iteration, the sampler updates the patient-specific activation
# probabilities, spike-and-slab parameters, observation probabilities, loading
# matrices, shared pattern matrix, and dependence parameters. Posterior samples
# retained after burn-in are returned for subsequent inference and uncertainty
# quantification.
#
# The function uses a "+1" storage convention: the supplied initial state is
# stored as the first sample, followed by the retained post-burn-in MCMC draws.
#
# Inputs:
#   X1       : A K × G binary diagnosis data matrix.
#   X2       : A K × G binary relapse data matrix.
#   W1_o     : Initial K × R binary diagnosis loading matrix.
#   W2_o     : Initial K × R binary relapse loading matrix.
#   H_o      : Initial R × G binary shared pattern matrix.
#   psi_o    : Initial length-G vector of spike-and-slab indicators.
#   n_iter   : Total number of Gibbs sampling iterations.
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
#   seed     : Optional random seed for reproducibility.
#
# Output:
#   A list containing posterior samples of:
#
#     samples_W1     : Diagnosis loading matrices.
#     samples_W2     : Relapse loading matrices.
#     samples_H      : Shared pattern matrices.
#     samples_alpha1 : Patient-specific activation probabilities for W1.
#     samples_beta   : Chromosome arm-specific inclusion probabilities.
#     samples_p1_11  : Diagnosis observation probabilities P(X1 = 1 |
#                       X1_tilde = 1).
#     samples_p1_10  : Diagnosis observation probabilities P(X1 = 1 |
#                       X1_tilde = 0).
#     samples_p2_11  : Relapse observation probabilities P(X2 = 1 |
#                       X2_tilde = 1).
#     samples_p2_10  : Relapse observation probabilities P(X2 = 1 |
#                       X2_tilde = 0).
#     samples_gammas : Posterior samples of gamma11 and gamma00.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
dependent_JBBMF <- function(X1, X2,W1_o, W2_o, H_o, psi_o,n_iter, burn_in,
                            a1, a2, b1, b2, c1, c2, d1, d2,u11, v11, u00, v00,seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # ----- Initialize state -----
  W1 <- W1_o
  W2 <- W2_o
  H  <- H_o
  psi_g <- psi_o
  K <- nrow(W1); R <- ncol(W1); G <- ncol(X1)
  # Initialize global gammas from priors (from counts)
  gamma11 <- rbeta(1, u11, v11)
  gamma00 <- rbeta(1, u00, v00)
  # ----- Preallocate (+1 scheme) -----
  num_samples <- n_iter - burn_in + 1
  if (num_samples <= 0) stop("burn_in must be < n_iter (or use the +1 scheme).")
  
  samples_W1      <- array(0L, dim = c(K, R, num_samples))
  samples_W2      <- array(0L, dim = c(K, R, num_samples))
  samples_H       <- array(0L, dim = c(R, G, num_samples))
  samples_alpha1  <- matrix(0, nrow = K, ncol = num_samples)
  samples_beta    <- matrix(0, nrow = G, ncol = num_samples)
  samples_p1_11   <- numeric(num_samples)
  samples_p1_10   <- numeric(num_samples)
  samples_p2_11   <- numeric(num_samples)
  samples_p2_10   <- numeric(num_samples)
  samples_gammas  <- matrix(NA_real_, nrow = num_samples, ncol = 2,
                            dimnames = list(NULL, c("gamma11","gamma00")))
  # ----- Save initial state as sample 1 -----
  # alpha, beta_g and p’s computed at the initial state
  alpha1_init <- update_alpha(W1, a1, a2)
  beta_g_init <- update_beta_g(H, psi_g, b1, b2, c1, c2)
  p1_init <- update_p(W1, H, X1)  # returns c(p11, p10)
  p2_init <- update_p(W2, H, X2)
  
  samples_W1[,,1]     <- W1
  samples_W2[,,1]     <- W2
  samples_H[,,1]      <- H
  samples_alpha1[,1]  <- alpha1_init
  samples_beta[,1]    <- beta_g_init
  samples_p1_11[1]    <- p1_init[1]; samples_p1_10[1] <- p1_init[2]
  samples_p2_11[1]    <- p2_init[1]; samples_p2_10[1] <- p2_init[2]
  samples_gammas[1, ] <- c(gamma11, gamma00)
  
  # ----- Progress bar -----
  pb <- .progress_bar_or_null(n_iter)
  
  # ----- MCMC iterations -----
  for (iter in 1:n_iter) {
    
    # (1) Update hyper/aux blocks not depending on the "new" W’s
    alpha1 <- update_alpha(W1, a1, a2)
    pi_s   <- update_pi(psi_g, d1, d2)
    betas  <- update_beta(H, b1, b2)  # if you use per-g Beta(prior) for H entries
    beta_g <- update_beta_g(H, psi_g, b1, b2, c1, c2)
    psi_g  <- update_psi(beta_g, b1, b2, c1, c2, pi_s)
    
    # (2) Update noise parameters p_{i,11}, p_{i,10}
    pvals1 <- update_p(W1, H, X1); p1_11 <- pvals1[1]; p1_10 <- pvals1[2]
    pvals2 <- update_p(W2, H, X2); p2_11 <- pvals2[1]; p2_10 <- pvals2[2]
    
    # (3) Update W1 using its prior alpha1 and view-1 likelihood
    W1 <- update_W(X1, W1, H, alpha = alpha1, p11 = p1_11, p10 = p1_10)
    
    # (4) Update W2 using **current global** gammas (do NOT resample here)
    W2 <- update_W2_Depend(W2, W1, H, X2, p11=p2_11, p10=p2_10,gamma11,gamma00, seed = NULL)
    
    # (5) Update H using both X1 and X2
    H <- update_H(X1, X2, W1, W2, H,p1_11, p1_10, p2_11, p2_10,beta = betas)
    
    # (6) Now resample gammas **once** using UPDATED W1, W2
    gams <- gammas_W1_W2(W1, W2, u11, v11, u00, v00)  # returns c(gamma11=., gamma00=.)
    gamma11 <- unname(gams["gamma11"])
    gamma00 <- unname(gams["gamma00"])
    
    if (!is.null(pb)) pb$tick()
    # (7) Store if beyond burn-in (using +1 scheme index)
    if (iter > burn_in) {
      idx <- iter - burn_in + 1  # with burn_in=0, iter=1 -> idx=2
      samples_W1[,,idx]     <- W1
      samples_W2[,,idx]     <- W2
      samples_H[,,idx]      <- H
      samples_alpha1[,idx]  <- alpha1
      samples_beta[,idx]    <- beta_g
      samples_p1_11[idx]    <- p1_11
      samples_p1_10[idx]    <- p1_10
      samples_p2_11[idx]    <- p2_11
      samples_p2_10[idx]    <- p2_10
      samples_gammas[idx, ] <- c(gamma11, gamma00)
    }
  }
  
  if (!is.null(pb)) pb$terminate()
  return(list(
    samples_W1    = samples_W1,
    samples_W2    = samples_W2,
    samples_H     = samples_H,
    samples_alpha1= samples_alpha1,
    samples_beta  = samples_beta,
    samples_p1_11 = samples_p1_11,
    samples_p1_10 = samples_p1_10,
    samples_p2_11 = samples_p2_11,
    samples_p2_10 = samples_p2_10,
    samples_gammas= samples_gammas
  ))
}