###############################################################################
# noise_bool.R
#
# Introduces random bit-flip noise into a binary matrix.
#
# This function randomly flips a specified proportion of entries in a binary
# matrix. Selected 0s are changed to 1s and selected 1s are changed to 0s,
# generating a noisy version of the original matrix. The function is primarily
# intended for simulation studies to assess the robustness of the JBBMF model
# under varying levels of observation noise.
#
# Inputs:
#   matrix      : A binary matrix with entries in {0,1}.
#   noise_level : Proportion of matrix entries to flip. Must lie between
#                 0 and 1.
#   seed        : Optional random seed for reproducibility.
#
# Output:
#   A binary matrix of the same dimensions as the input with the specified
#   proportion of entries randomly flipped.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################
# Function to introduce noise to a boolean matrix
noise_bool <- function(matrix, noise_level, seed=NULL) {
  # Check if noise_level is between 0 and 1
  if (noise_level < 0 || noise_level > 1) {
    stop("noise_level must be between 0 and 1")
  }
  # Determine the number of elements to flip
  num_elements <- length(matrix)
  num_flips <- ceiling(noise_level * num_elements)
  # Randomly select elements to flip
  # Setting the seed for reproducibility
  # Optionally set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  flip_indices <- sample(seq_len(num_elements), num_flips, replace = FALSE)
  # Flip the selected elements
  matrix[flip_indices] <- 1 - matrix[flip_indices]
  return(matrix)
}
