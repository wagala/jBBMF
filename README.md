# Joint Bayesian Boolean Matrix Factorization (JBBMF)

## Requirements

The implementation is written entirely in **R**.

The following packages are required:

```r
library(Matrix)
library(rBMF)
library(progress)
library(viridis)
library(pbmcapply)
```

These packages provide support for

- **Matrix** – efficient sparse and dense matrix operations.
- **rBMF** – Boolean matrix factorization utilities.
- **progress** – progress bars for long-running computations.
- **viridis** – perceptually uniform color palettes for data visualization.
- **pbmcapply** – parallel execution with progress bars for multi-chain MCMC sampling.
## Overview

This repository contains the R implementation of the **Joint Bayesian Boolean Matrix Factorization (JBBMF)** model for jointly analyzing paired binary datasets with a shared latent structure.

The model was developed to analyze paired genomic datasets (e.g., diagnosis and relapse samples in multiple myeloma), where each dataset is represented by its own loading matrix while sharing a common latent pattern matrix.

Unlike standard joint Boolean matrix factorization approaches that assume the two loading matrices are independent, this implementation models the dependence between them through conditional transition probabilities.

---

## Model

Given two binary data matrices

- **Diagnosis:** \(X_1 \in \{0,1\}^{K\times G}\)
- **Relapse:** \(X_2 \in \{0,1\}^{K\times G}\)

the latent Boolean representations are

\[
\tilde X_1 = W_1 \circ H,
\]

\[
\tilde X_2 = W_2 \circ H,
\]

where

- \(W_1\) is the diagnosis loading matrix,
- \(W_2\) is the relapse loading matrix,
- \(H\) is the shared binary pattern matrix,
- \(\circ\) denotes the Boolean matrix product.

Dependence between the loading matrices is introduced through

\[
\gamma_{11}=P(W_2=1\mid W_1=1),
\]

\[
\gamma_{00}=P(W_2=0\mid W_1=0),
\]

allowing latent factors to persist or disappear between diagnosis and relapse.

Posterior inference is performed using a Gibbs sampler.

---

## Repository Structure

```
.
├── R/
│   ├── boolean_product.R
│   ├── computeLoglikelW.R
│   ├── computeLoglikelihoodTotal.R
│   ├── counts_for_W1_W2.R
│   ├── dependent_JBBMF.R
│   ├── gammas_W1_W2.R
│   ├── log_W2_given_W1.R
│   ├── logpost_W2_W1.R
│   ├── mult_dep_JBBMF.R
│   ├── update_alpha.R
│   ├── update_beta.R
│   ├── update_beta_g.R
│   ├── update_H.R
│   ├── update_p.R
│   ├── update_pi.R
│   ├── update_psi.R
│   ├── update_W.R
│   ├── update_W2_Depend.R
│   └── ...
├── data/
├── simulations/
├── results/
└── README.md
```

---

## Main Components

The implementation includes routines for

- Boolean matrix multiplication
- Gibbs sampling updates for \(W_1\), \(W_2\), and \(H\)
- Estimation of observation error probabilities
- Spike-and-slab priors for chromosome arm inclusion probabilities
- Estimation of dependence parameters \(\gamma_{11}\) and \(\gamma_{00}\)
- Multi-chain MCMC sampling
- Posterior summaries and uncertainty quantification

---

## Methodology

The implementation follows the Bayesian model described in the accompanying manuscript:

> **Joint Bayesian Boolean Matrix Factorization for Paired Binary Data with Dependent Latent Structures**

Posterior inference is performed using Gibbs sampling with conjugate updates whenever available.

---

## Requirements

The implementation is written entirely in **R**.

Required packages include

- Matrix
- parallel
- pbmcapply

Additional visualization packages may be used for posterior summaries and figures.

---

## Status

This repository is under active development.

Current work includes

- methodological development,
- simulation studies,
- real-data analysis,
- manuscript preparation.

---

## Author

**Adolphus Wagala**

Department of Biostatistics

Dana-Farber Cancer Institute / Harvard T.H. Chan School of Public Health

---

## License

This code is provided for research purposes.

A formal software license will be added upon public release.