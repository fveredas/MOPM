# MOPM.R (METHIONINE OXIDATION PREDICTIVE MODEL)
# 
# Copyright (c) 2017 Juan Carlos Aledo & Francisco J. Veredas
# -------------------------------------------------------------------------------
# -------------------- MIT License ----------------------------------------------
# -------------------------------------------------------------------------------
# Permision is hereby granted, free of charge, to any person obtaining a copy of 
# this software and associated documentation files (the "Software", to deal in the
# Software without restriction), including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
# Software, and to permitto whom the Software is furnished to do so, subjecto to
# the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all
# copies or substanctial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------
# --------------- Brief summary of what this sofware does --------------------------
# ----------------------------------------------------------------------------------
# The current program is an implimentation of the predictive model described in 
# Neural Computing & Applications... It takes as input a PDB ID (the 4-character unique 
# identifier of every entry in the Protein Data Bank) and return a dataframe 
# with as many rows as methionine residues are found in the protein being anlyzed.
# For each methionine, this dataframe makes a prediction with respect to its 
# oxidation status ("Yes" if the residue is predicted to be oxidezed, and "No" otherwise). 
# Each prediction is accompanied by its probability. In addition, the value of the
# extrated features is also provided. The dataframe is saved as a Rda file in the user 
# computer with the name "pdbID_MOPM.Rda".
# (Note: this script has been optimised to be run at Unix computers. It also works at
# MS Windows computers, although for a clean execution "system" in the .R scripts
# should be changed for the corresponding commented options for Windows).
# ------------------------------------------------------------------------------------

# ----------------------------- Dependencies 
# Make sure you have installed the following R packages:
# install.packages("bio3d")
# install.packages("httr")
# install.packages("Rpdb")
# install.packages("RCurl")
# install.packages("bitops")
# install.packages("caret")
# install.packages("randomForest")
# install.packages("ipred")
# In addition, you will need to install MUSCLE and DSSP. To this end,
# you may wish to follow the indications given at:
# http://thegrantlab.org/bio3d/tutorials/installing-bio3d

## ---------------- Working Directory
# Set the sorking directory
# In the same working directory where the current script is located, there must be found also 
# the following files:
# MOPM.R <-- main R script
# Methionine_Oxidation_Predictive_Model.R
# ToolBox.R
# dif.R
# pdb2acc.R
# Phosphorylation_site_dataset
# btaurus_orth.rda
# drerio_orth.rda
# ggallus_orth.rda
# ggorilla_orth.rda
# prdVars_log_SBF.RData
# prdVars_log_ext_SBF.RData
# preProc_final.Rdata
# ptroglodytes_orth.rda
# rf_final.Rda
# rnorvegicus_orth.rda
# xtropicalis_orth.rda


rm(list=ls())

source("Methionine_Oxidation_Predictive_Model.R")

MOPM_results <- predict.MetO(pdbID = "3CWM",	# Protein's pdb identifier
                chain = "A",					# Protein's chain to be analysed
                email = "user@server", 			# Your email
                cut = 0.5040)					# Probability cutoff for "oxidised" prediction