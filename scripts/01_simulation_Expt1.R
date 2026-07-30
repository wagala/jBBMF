###############################################################################
# simulation_scenario1.R
#
# Loads all JBBMF functions and generates the latent matrices and binary data
# for Scenario 1: relapse gains additional latent factors.
#
# Run this script from the project root directory.
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################

#===============================================================================
# Load all functions
#===============================================================================

source("R/source_all.R")


#===============================================================================
# General simulation settings
#===============================================================================

set.seed(42)

K <- 62

# Chromosome arms: 1p, 1q, ..., 22p, 22q
arms <- unlist(lapply(1:22, function(i) paste0(i, c("p", "q"))))

G <- length(arms)  # 44

rownames_template <- paste0("X", 1:K)
colnames_template <- arms


#===============================================================================
# Scenario 1: Relapse gains extra factors
#===============================================================================

R_shared  <- 2
R_relapse <- 2
R <- R_shared + R_relapse


#-------------------------------------------------------------------------------
# Construct the shared pattern matrix H
#-------------------------------------------------------------------------------

# The first R_shared rows represent shared block factors.
# The remaining rows represent relapse-specific dispersed factors.

H <- matrix(0L, R, G)

cols_shared <- split(
  1:G,
  cut(1:G, R_shared, labels = FALSE)
)

for (r in 1:R_shared) {
  H[r, cols_shared[[r]]] <- 1L
}


#-------------------------------------------------------------------------------
# Alternative specification: Experiment 1b
#-------------------------------------------------------------------------------

# cols1 <- c(3, 6, 9, 12, 14, 17)
# cols2 <- c(23, 31, 32, 33, 34, 35, 38, 40, 41, 42, 43, 44)
# H[4, cols1] <- 1L
# H[4, cols2] <- 1L


#-------------------------------------------------------------------------------
# Experiment I: relapse-specific dispersed factors
#-------------------------------------------------------------------------------

set.seed(101)

for (r in (R_shared + 1):R) {
  H[r, sample(1:G, 10)] <- 1L
}


#-------------------------------------------------------------------------------
# Construct the diagnosis loading matrix W1
#-------------------------------------------------------------------------------

# Patients are exposed only to the shared factors at diagnosis.

W1 <- matrix(0L, K, R)

for (r in 1:R_shared) {
  W1[sample(1:K, 18), r] <- 1L
}


#-------------------------------------------------------------------------------
# Construct the relapse loading matrix W2
#-------------------------------------------------------------------------------

# Relapse retains the diagnosis exposures and adds relapse-specific exposures.

W2 <- W1

for (r in (R_shared + 1):R) {
  W2[sample(1:K, 12), r] <- 1L
}


#-------------------------------------------------------------------------------
# Generate the latent Boolean data matrices
#-------------------------------------------------------------------------------

X1 <- boolean_product(W1, H)
X2 <- boolean_product(W2, H)

# Alternative Boolean product:
#
# UU1 <- as.matrix(W1 %&% H) * 1
# UU2 <- as.matrix(W2 %&% H) * 1


#-------------------------------------------------------------------------------
# Assign row and column names
#-------------------------------------------------------------------------------

rownames(X1) <- rownames(X2) <- rownames_template
colnames(X1) <- colnames(X2) <- colnames_template

rownames(W1) <- rownames_template
rownames(W2) <- rownames_template

colnames(H) <- colnames_template
#===============================================================================
# Reshape the true matrices for visualization
#===============================================================================

W1_melt <- reshape2::melt(W1)
W1_melt$Var2 <- as.factor(W1_melt$Var2)

W2_melt <- reshape2::melt(W2)
W2_melt$Var2 <- as.factor(W2_melt$Var2)

H_melt <- reshape2::melt(H)
H_melt$Var1 <- as.factor(H_melt$Var1)

X1_melt <- reshape2::melt(X1)
X2_melt <- reshape2::melt(X2)


#===============================================================================
# Introduce observation noise
#===============================================================================

noise_level <- 0.15

X1_noisy <- noise_bool(
  matrix      = X1,
  noise_level = noise_level,
  seed        = 123
)

X2_noisy <- noise_bool(
  matrix      = X2,
  noise_level = noise_level,
  seed        = 123
)


#===============================================================================
# Reshape the noisy data matrices for visualization
#===============================================================================

X1_melt_noisy <- reshape2::melt(X1_noisy)
X2_melt_noisy <- reshape2::melt(X2_noisy)
#===============================================================================
# Prepare melted matrices for plotting
#===============================================================================

# Preserve the desired factor ordering.
W1_melt$Var2 <- factor(W1_melt$Var2, levels = 1:R)
W2_melt$Var2 <- factor(W2_melt$Var2, levels = 1:R)

# Reverse the factor levels so that Factor 1 appears at the top.
H_melt$Var1 <- factor(H_melt$Var1, levels = R:1)

# Preserve chromosome arm ordering.
X1_melt$Var2 <- factor(X1_melt$Var2, levels = arms)
X2_melt$Var2 <- factor(X2_melt$Var2, levels = arms)

X1_melt_noisy$Var2 <- factor(
  X1_melt_noisy$Var2,
  levels = arms
)

X2_melt_noisy$Var2 <- factor(
  X2_melt_noisy$Var2,
  levels = arms
)

# Convert binary values to factors for discrete legends.
W1_melt$value <- factor(W1_melt$value, levels = c(0, 1))
W2_melt$value <- factor(W2_melt$value, levels = c(0, 1))
H_melt$value  <- factor(H_melt$value,  levels = c(0, 1))

X1_melt$value <- factor(X1_melt$value, levels = c(0, 1))
X2_melt$value <- factor(X2_melt$value, levels = c(0, 1))

X1_melt_noisy$value <- factor(
  X1_melt_noisy$value,
  levels = c(0, 1)
)

X2_melt_noisy$value <- factor(
  X2_melt_noisy$value,
  levels = c(0, 1)
)


#===============================================================================
# Common theme for binary heatmaps
#===============================================================================

binary_heatmap_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    
    # Place the discrete legend on the right.
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    # Reduce axis font sizes.
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 5),
    
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    
    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 9,
      hjust = 0.5
    )
  )


#===============================================================================
# Diagnosis loading matrix W1
#===============================================================================

plot_W1 <- ggplot(
  W1_melt,
  aes(x = Var2, y = Var1, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.15
  ) +
  scale_fill_viridis_d(
    option = "C",
    limits = c("0", "1"),
    breaks = c("0", "1"),
    labels = c("0", "1"),
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = paste0("X", nrow(W1):1),
    drop = FALSE,
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  coord_fixed() +
  labs(
    title = expression("Diagnosis loading matrix " * bold(W)[1]),
    x = "Factor",
    y = "Patient"
  ) +
  binary_heatmap_theme

plot_W1


#===============================================================================
# Relapse loading matrix W2
#===============================================================================

plot_W2 <- ggplot(
  W2_melt,
  aes(x = Var2, y = Var1, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.15
  ) +
  scale_fill_viridis_d(
    option = "C",
    limits = c("0", "1"),
    breaks = c("0", "1"),
    labels = c("0", "1"),
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = paste0("X", nrow(W1):1),
    drop = FALSE,
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  coord_fixed() +
  labs(
    title = expression("Relapse loading matrix " * bold(W)[2]),
    x = "Factor",
    y = "Patient"
  ) +
  binary_heatmap_theme

plot_W2


#===============================================================================
# Shared pattern matrix H
#===============================================================================

plot_H <- ggplot(
  H_melt,
  aes(x = Var2, y = Var1, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.15
  ) +
  scale_fill_viridis_d(
    option = "C",
    limits = c("0", "1"),
    breaks = c("0", "1"),
    labels = c("0", "1"),
    drop = FALSE
  ) +
  scale_x_discrete(
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    labels = paste0( R:1),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  labs(
    title = expression("Shared pattern matrix " * bold(H)),
    x = "Chromosome arm",
    y = "Factor"
  ) +
  binary_heatmap_theme +
  theme(
    axis.text.x = element_text(
      size = 6,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    axis.text.y = element_text(size = 7)
  )

plot_H


#===============================================================================
# True diagnosis matrix X1
#===============================================================================

plot_X1 <- ggplot(
  X1_melt,
  aes(x = Var2, y = Var1, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.10
  ) +
  scale_fill_viridis_d(
    option = "C",
    limits = c("0", "1"),
    breaks = c("0", "1"),
    labels = c("0", "1"),
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = paste0("X", nrow(W1):1),
    drop = FALSE,
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    limits = arms,
    expand = c(0, 0)
  ) +
  coord_fixed() +
  labs(
    title = expression("True diagnosis matrix " * bold(X)[1]),
    x = "Chromosome arm",
    y = "Patient"
  ) +
  binary_heatmap_theme +
  theme(
    axis.text.x = element_text(
      size = 5,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    axis.text.y = element_text(size = 4)
  )

plot_X1


#===============================================================================
# True relapse matrix X2
#===============================================================================

plot_X2 <- ggplot(
  X2_melt,
  aes(x = Var2, y = Var1, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.10
  ) +
  scale_fill_viridis_d(
    option = "C",
    limits = c("0", "1"),
    breaks = c("0", "1"),
    labels = c("0", "1"),
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = paste0("X", nrow(W1):1),
    drop = FALSE,
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    limits = arms,
    expand = c(0, 0)
  ) +
  coord_fixed() +
  labs(
    title = expression("True relapse matrix " * bold(X)[2]),
    x = "Chromosome arm",
    y = "Patient"
  ) +
  binary_heatmap_theme +
  theme(
    axis.text.x = element_text(
      size = 5,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    axis.text.y = element_text(size = 4)
  )

plot_X2


#===============================================================================
# Noisy diagnosis matrix X1
#===============================================================================

plot_X1_noisy <- ggplot(
  X1_melt_noisy,
  aes(x = Var2, y = Var1, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.10
  ) +
  scale_fill_viridis_d(
    option = "C",
    limits = c("0", "1"),
    breaks = c("0", "1"),
    labels = c("0", "1"),
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = paste0("X", nrow(W1):1),
    drop = FALSE,
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    limits = arms,
    expand = c(0, 0)
  ) +
  coord_fixed() +
  labs(
    title = expression("Noisy diagnosis matrix " * bold(X)[1]),
    subtitle = "Noise level = 0.15",
    x = "Chromosome arm",
    y = "Patient"
  ) +
  binary_heatmap_theme +
  theme(
    axis.text.x = element_text(
      size = 5,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    axis.text.y = element_text(size = 4)
  )

plot_X1_noisy


#===============================================================================
# Noisy relapse matrix X2
#===============================================================================

plot_X2_noisy <- ggplot(
  X2_melt_noisy,
  aes(x = Var2, y = Var1, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.10
  ) +
  scale_fill_viridis_d(
    option = "C",
    limits = c("0", "1"),
    breaks = c("0", "1"),
    labels = c("0", "1"),
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = paste0("X", nrow(W1):1),
    drop = FALSE,
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    limits = arms,
    expand = c(0, 0)
  ) +
  coord_fixed() +
  labs(
    title = expression("Noisy relapse matrix " * bold(X)[2]),
    subtitle = "Noise level = 0.15",
    x = "Chromosome arm",
    y = "Patient"
  ) +
  binary_heatmap_theme +
  theme(
    axis.text.x = element_text(
      size = 5,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 4),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm")
  )

plot_X2_noisy
#===============================================================================
# Prepare observed data
#===============================================================================

# Use the noisy diagnosis and relapse matrices as the observed datasets.
X11 <- X1_noisy
X22 <- X2_noisy

# Convert to logical matrices required by ASSO.
Xb1 <- X11 == 1
Xb2 <- X22 == 1

# Number of chromosome arms.
G <- ncol(X11)

# Number of latent factors.
R <- 4


#===============================================================================
# Obtain initial values using ASSO
#===============================================================================

# ASSO provides deterministic initial estimates of the loading and
# pattern matrices to initialize the Gibbs sampler.

res_asso1 <- Asso_approximate(
  Xb1,
  R,
  list(
    threshold = 0.5,
    penalty_overcovered = 1,
    bonus_covered = 1,
    verbose = 0
  )
)

res_asso2 <- Asso_approximate(
  Xb2,
  R,
  list(
    threshold = 0.5,
    penalty_overcovered = 1,
    bonus_covered = 1,
    verbose = 0
  )
)


#===============================================================================
# Initial values
#===============================================================================

# Diagnosis loading matrix.
W1_o <- as.matrix(res_asso1$O) * 1

# Relapse loading matrix.
W2_o <- as.matrix(res_asso2$O) * 1

# Shared pattern matrix.
H_o <- as.matrix(res_asso1$B) * 1

# Individual ASSO pattern matrices retained for comparison.
H1 <- as.matrix(res_asso1$B) * 1
H2 <- as.matrix(res_asso2$B) * 1

# Initial spike-and-slab indicators.
psi_o <- rbinom(
  G,
  size = 1,
  prob = 0.5
)


#===============================================================================
# Prior hyperparameters
#===============================================================================

# Beta prior for chromosome-arm inclusion probabilities.
c2 <- 1
d2 <- 1

# Beta priors governing dependence between W1 and W2.
uu11 <- 1
vv11 <- 1

uu00 <- 1
vv00 <- 1


#===============================================================================
# MCMC settings
#===============================================================================

burn_in <- 0

n_chains <- 4

n <- 1000


#===============================================================================
# Fit the dependent JBBMF model
#===============================================================================
# Run one chain without parallel processing
test_chain <- dependent_JBBMF(
  X1 = X11,
  X2 = X22,
  W1_o = W1_o,
  W2_o = W2_o,
  H_o = H_o,
  psi_o = psi_o,
  n_iter = 10,
  burn_in = 0,
  a1 = 1,
  a2 = 1,
  b1 = 1,
  b2 = 1,
  c1 = 1,
  c2 = c2,
  d1 = 1,
  d2 = d2,
  u11 = uu11,
  v11 = vv11,
  u00 = uu00,
  v00 = vv00,
  seed = 1357
)
### Run 4 chains in parrallel
sim_ExptI <- mult_dep_JBBMF(
  X1 = X11,
  X2 = X22,
  W1_o = W1_o,
  W2_o = W2_o,
  H_o = H_o,
  psi_o = psi_o,
  n_iter = n,
  burn_in = burn_in,
  a1 = 1,
  a2 = 1,
  b1 = 1,
  b2 = 1,
  c1 = 1,
  c2 = c2,
  d1 = 1,
  d2 = d2,
  u11 = uu11,
  v11 = vv11,
  u00 = uu00,
  v00 = vv00,
  seed = 1357,
  n_chains = n_chains
)
