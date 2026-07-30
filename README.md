# Joint Bayesian Boolean Matrix Factorization (JBBMF)

## Overview

This repository contains the R implementation of the **Joint Bayesian Boolean Matrix Factorization (JBBMF)** model for jointly analyzing paired binary datasets with a shared latent structure.

The methodology was developed for paired genomic studies, such as diagnosis and relapse samples in multiple myeloma, where each dataset is represented by its own latent loading matrix while sharing a common binary pattern matrix. Unlike existing joint Boolean matrix factorization approaches that assume independent loading matrices, JBBMF explicitly models their dependence through conditional transition probabilities.

---

## Model

Let

- $\mathbf{X}_1,\mathbf{X}_2\in\{0,1\}^{K\times G}$ denote paired binary data matrices,
- $\mathbf{W}_1,\mathbf{W}_2\in\{0,1\}^{K\times R}$ denote the corresponding latent loading matrices,
- $\mathbf{H}\in\{0,1\}^{R\times G}$ denote the shared binary pattern matrix.

The latent Boolean representations are

$$
\tilde{\mathbf{X}}_1=\mathbf{W}_1\circ\mathbf{H},
\qquad
\tilde{\mathbf{X}}_2=\mathbf{W}_2\circ\mathbf{H},
$$

where $\circ$ denotes the Boolean matrix product.

Dependence between diagnosis and relapse is modeled through

$$
\gamma_{11}=P\!\left(W_{2,kr}=1 \mid W_{1,kr}=1\right),
$$

and

$$
\gamma_{00}=P\!\left(W_{2,kr}=0\mid W_{1,kr}=0 \right),
$$

allowing latent factors to persist or disappear between paired samples.

Posterior inference is performed using a Gibbs sampler.

---

## Repository Structure

```text
.
├── R/
├── scripts/
├── data/
├── results/
├── README.md
└── jBBMF.Rproj
```

- **R/** – Core functions implementing the Bayesian model and Gibbs sampler.
- **scripts/** – Reproducible scripts for simulations and data analyses.
- **data/** – Input datasets.
- **results/** – Generated figures, tables, and posterior summaries.

---

## Main Components

The implementation includes routines for

- Boolean matrix multiplication
- Gibbs sampling for $\mathbf{W}_1$, $\mathbf{W}_2$, and $\mathbf{H}$
- Estimation of observation error probabilities
- Spike-and-slab priors for feature inclusion
- Estimation of the dependence parameters $\gamma_{11}$ and $\gamma_{00}$
- Multi-chain MCMC sampling
- Posterior summaries and uncertainty quantification

---

## Requirements

The implementation is written entirely in **R**.

Required packages:

```r
library(Matrix)
library(rBMF)
library(progress)
library(viridis)
library(pbmcapply)
```

These packages provide

- **Matrix** – efficient sparse and dense matrix operations.
- **rBMF** – Boolean matrix factorization utilities.
- **progress** – progress bars for long-running computations.
- **viridis** – perceptually uniform colour palettes for visualization.
- **pbmcapply** – parallel execution with progress bars for multi-chain MCMC sampling.

---

## Running the Code

Clone the repository and source all functions

```r
source("R/source_all.R")
```

Simulation studies and analyses can then be reproduced by running the scripts in the `scripts/` directory.

---

## Status

This software is under active development and accompanies the manuscript

> **Joint Bayesian Boolean Matrix Factorization for Paired Binary Data with Dependent Latent Structures**

Future releases will include additional simulation studies, real-data analyses, and an installable R package.

---

## Author

**Adolphus Wagala**  
Dept. of Data Science and Dept. of Biostatistics  
Dana-Farber Cancer Institute / Harvard T.H. Chan School of Public Health