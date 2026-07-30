###############################################################################
# runJBBMF.R
#
# Runs the Gibbs sampler for the sparse joint Bayesian Boolean matrix
# factorization model.
#
# The function jointly analyzes paired binary diagnosis and relapse datasets,
# X1 and X2. The two datasets have separate patient loading matrices, W1 and
# W2, but share the pattern matrix H:
#
#   X1_tilde = W1 ∘ H,
#   X2_tilde = W2 ∘ H,
#
# where ∘ denotes the Boolean matrix product.
#
# At each Gibbs iteration, the sampler updates:
#
#   1. Patient-specific activation probabilities for W1 and W2;
#   2. The spike-and-slab mixing probability;
#   3. Chromosome arm-specific inclusion probabilities;
#   4. Spike-and-slab allocation indicators;
#   5. Observation model parameters for X1 and X2;
#   6. The diagnosis loading matrix W1;
#   7. The relapse loading matrix W2; and
#   8. The shared pattern matrix H.
#
# Element-wise Gibbs updates are used for W1, W2, and H. The function stores
# their retained samples together with samples of the patient-specific
# activation probabilities, chromosome arm-specific inclusion probabilities,
# and observation model parameters.
#
# Inputs:
#   X1       : A K × G binary diagnosis data matrix.
#   X2       : A K × G binary relapse data matrix.
#   W1_o     : An initial K × R binary diagnosis loading matrix.
#   W2_o     : An initial K × R binary relapse loading matrix.
#   H_o      : An initial R × G binary shared pattern matrix.
#   psi_o    : A binary vector of length G containing the initial spike-and-
#              slab allocation indicators for the chromosome arms.
#   n_iter   : Total number of Gibbs sampling iterations.
#   burn_in  : Number of initial Gibbs iterations discarded as burn-in.
#   a1, a2   : Shape parameters of the Beta prior for the patient-specific
#              activation probabilities of W1 and W2.
#   b1, b2   : Shape parameters of the slab Beta distribution for the
#              chromosome arm-specific inclusion probabilities.
#   c1, c2   : Shape parameters of the spike Beta distribution for the
#              chromosome arm-specific inclusion probabilities.
#   d1, d2   : Shape parameters of the Beta prior for the spike-and-slab
#              mixing probability.
#
# Output:
#   A named list containing:
#
#   samples_p1_10    : Samples of P(X1 = 1 | X1_tilde = 0).
#   samples_p1_11    : Samples of P(X1 = 1 | X1_tilde = 1).
#   samples_p2_10    : Samples of P(X2 = 1 | X2_tilde = 0).
#   samples_p2_11    : Samples of P(X2 = 1 | X2_tilde = 1).
#   samples_alpha1   : Samples of the diagnosis patient-specific activation
#                      probabilities.
#   samples_alpha2   : Samples of the relapse patient-specific activation
#                      probabilities.
#   samples_beta     : Samples of the chromosome arm-specific inclusion
#                      probabilities under the spike-and-slab prior.
#   samples_W1       : Retained samples of the diagnosis loading matrix.
#   samples_W2       : Retained samples of the relapse loading matrix.
#   samples_H        : Retained samples of the shared pattern matrix.
#
# Notes:
#   The first position in each returned object contains an initialized value.
#   Subsequent positions contain draws obtained after the burn-in period.
#
#   The function requires the progress package because it uses
#   progress_bar$new() to display sampling progress.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
#-------------------------Gibbs Sampler for Sparse BBMF-------------------------
##Code that stores the initial values
runJBBMF<- function(X1,X2, W1_o,W2_o, H_o, psi_o, n_iter, burn_in,
                    a1, a2, b1, b2, c1, c2, d1, d2) {
  # Initialize Ws and H with the initial values provided
  W1 <- W1_o
  W2 <- W2_o
  H <- H_o
  # Dimensions for W and H
  K <- nrow(W1_o)
  R <- ncol(W1_o)
  G <- ncol(X1)
  psi_g<-psi_o
  # Pre-allocate space for samples, including space for initial values
  num_samples <- n_iter - burn_in + 1
  samples_W1 <- array(0, dim = c(K, R, num_samples))
  samples_W2 <- array(0, dim = c(K, R, num_samples))
  samples_H <- array(0, dim = c(R, G, num_samples))
  samples_alpha1 <- matrix(0, nrow = K, ncol = num_samples)
  samples_alpha2 <- matrix(0, nrow = K, ncol = num_samples)
  samples_beta <- matrix(0, nrow = G, ncol = num_samples)
  samples_p1_11 <- numeric(num_samples)
  samples_p1_10 <- numeric(num_samples)
  samples_p2_11 <- numeric(num_samples)
  samples_p2_10 <- numeric(num_samples)
  
  # Store initial values
  samples_W1[,,1] <- W1
  samples_W2[,,1] <- W2
  samples_H[,,1] <- H
  samples_alpha1[,1] <- update_alpha(W1, a1, a2)  
  samples_alpha2[,1] <- update_alpha(W2, a1, a2)  
  samples_beta[,1] <- update_beta_g(H, psi_o, b1, b2, c1, c2) 
  pvals1 <- update_p(W1, H, X1)
  pvals2 <- update_p(W2, H, X2)
  samples_p1_11[1] <- pvals1[1]
  samples_p1_10[1] <- pvals1[2]
  samples_p2_11[1] <- pvals2[1]
  samples_p2_10[1] <- pvals2[2]
  
  # Initialize the progress bar
  pb_chain <- progress_bar$new(format = "working [:bar] :percent in :elapsed",
                               total = n_iter, clear = FALSE, width = 100)
  
  for (iter in 1:n_iter) {
    # Perform updates for W, H, and other parameters
    alpha1 <- update_alpha(W1, a1, a2)
    alpha2 <- update_alpha(W2, a1, a2)
    pi_s <- update_pi(psi_g, d1, d2)
    betas<-update_beta(H,b1,b2)
    beta_g <- update_beta_g(H, psi_g, b1, b2, c1, c2)
    psi_g <- update_psi(beta_g, b1, b2, c1, c2, pi_s)
    pvals1 <- update_p(W1, H, X1)
    pvals2 <- update_p(W2, H, X2)
    p1_11 <- pvals1[1]
    p1_10 <- pvals1[2]
    p2_11 <- pvals2[1]
    p2_10 <- pvals2[2]
    
    # Update W and H
    W1<-update_W(X1,W1, H, alpha=alpha1, p11=p1_11, p10=p1_10)
    W2<-update_W(X2,W2, H, alpha=alpha2, p11=p2_11, p10=p2_10)
    H<-update_H(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta=betas) 
    # W1<-update_W_Block_rows(X1,W1,H, alpha=alpha1, p11=p1_11, p10=p1_10)
    # W2<-update_W_Block_rows(X2,W2,H, alpha=alpha2, p11=p2_11, p10=p2_10)
    #update_H_Block_columns(X1, X2, W1, W2, H,p1_11, p1_10,p2_11, p2_10, beta=beta_g)
    # Update the progress bar
    pb_chain$tick()
    # Store the results
    if (iter > burn_in) {
      idx <- iter - burn_in + 1
      samples_W1[,,idx] <- W1
      samples_W2[,,idx] <- W2
      samples_H[,,idx] <- H
      samples_alpha1[,idx] <- alpha1
      samples_alpha2[,idx] <- alpha2
      samples_beta[,idx] <- beta_g
      samples_p1_11[idx] <- p1_11
      samples_p1_10[idx] <- p1_10
      samples_p2_11[idx] <- p2_11
      samples_p2_10[idx] <- p2_10
    }
  }
  
  # Terminate the progress bar
  pb_chain$terminate()
  
  #Return a List of Arrays
  return(list(samples_p1_10=samples_p1_10,samples_p1_11=samples_p1_11,
              samples_p2_10=samples_p2_10,samples_p2_11=samples_p2_11,
              samples_alpha1=samples_alpha1,samples_alpha2=samples_alpha2,
              samples_beta=samples_beta,samples_W1=samples_W1,
              samples_W2=samples_W2,samples_H=samples_H
  ))  
}