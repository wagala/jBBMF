###############################################################################
# source_all.R
#
# Loads the required libraries and sources all functions used in the JBBMF
# analysis.
#
# Run from the project root using:
#
#   source("R/source_all.R")
#
# Author: Adolphus Wagala
# Created: 2026-07-30
###############################################################################


#===============================================================================
# Required libraries
#===============================================================================

library(Matrix)
library(rBMF)
library(progress)
library(viridis)
library(pbmcapply)
library(parallel)
#===============================================================================
# Simulation utilities
#===============================================================================

source("R/noise_bool.R")

#===============================================================================
# Boolean operations and numerical utilities
#===============================================================================

source("R/boolean_product.R")
source("R/sumLog.R")
#===============================================================================
# Likelihood functions
#===============================================================================
source("R/computeLoglikelW.R")

source("R/computeLoglikelihoodTotal.R")
#===============================================================================
# Prior and posterior calculations
#===============================================================================
source("R/compute_log_priorW.R")
source("R/compute_log_postW.R")

source("R/compute_log_priorH.R")
source("R/compute_log_posteriorH.R")
#===============================================================================
# Parameter updates
#===============================================================================
source("R/update_alpha.R")
source("R/update_beta.R")
source("R/update_beta_g.R")
source("R/update_p.R")
source("R/update_pi.R")
source("R/update_psi.R")
#===============================================================================
# Loading matrix updates
#===============================================================================
source("R/update_W.R")
source("R/update_W_Block_rows.R")
source("R/update_W2_Depend.R")
#===============================================================================
# Shared pattern matrix updates
#===============================================================================
source("R/update_H.R")
source("R/update_H_Block_columns.R")
#===============================================================================
# Dependence between W1 and W2
#===============================================================================
source("R/counts_for_W1_W2.R")
source("R/log_W2_given_W1.R")
source("R/gammas_W1_W2.R")
source("R/logpost_W2_W1.R")
#===============================================================================
# Independent JBBMF samplers
# Retained temporarily for comparison and validation.
#===============================================================================
source("R/runJBBMF.R")
source("R/run_multiple_JBBMF.R")
#===============================================================================
# Dependent JBBMF samplers
#===============================================================================
source("R/progress_bar_or_null.R")
source("R/dependent_JBBMF.R")
source("R/mult_dep_JBBMF.R")
source("R/dependent_JBBMF.R")
source("R/mult_dep_JBBMF.R")
#===============================================================================