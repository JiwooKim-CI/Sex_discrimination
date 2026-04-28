# Sex_discrimination

This repository contains code for analyzing sex discrimination using instrumental inequalities and partial-identification bounds. The project combines hypothesis testing for instrumental inequality violations with bound-based estimation of direct effects under an instrumental-variable framework.

## Overview

The repository includes two complementary analysis workflows:

1. **Instrumental inequality evidence analysis in R**
   - Computes point estimates of inequality evidence
   - Constructs one-sided Wald tests
   - Uses stratified bootstrap procedures to quantify uncertainty

2. **Bound-based analysis under IV assumptions**
   - Computes sharp bounds for direct-effect estimands such as
     - CDE(0)
     - CDE(1)
     - NDE(0)
     - NDE(1)
   - Implements both:
     - an R version using linear programming (`lpSolve`)
     - a Python version using `autobounds`

The code is designed around two empirical case studies:
- **Case 1:** promotion by sex and job type
- **Case 2:** regular vs. non-regular employment by sex and education

## Repository contents

- `sexdiscrimination_exclusion.R`  
  Main R script for instrumental inequality evidence analysis.  
  This script:
  - constructs the two datasets (Case 1 and Case 2),
  - computes the instrumental inequality point estimate,
  - runs Wald tests,
  - runs bootstrap resampling.

- `IV_bound.R`  
  R script for computing sharp bounds on direct effects under an IV-based potential outcomes model using linear programming.  
  It:
  - solves lower and upper bounds for different types of direct effects.

- `autobound.py`  
  Python script that uses the `autobounds` package to compute:
  - direct-effect bounds under the exclusion-restriction-violated model.

- `session_info.txt`  
  R session information for reproducibility, including R version and package versions.

## Data structure

The analyses use aggregated cell counts rather than individual-level raw data.

### Case 1
- `Z`: sex (`0 = female`, `1 = male`)
- `X`: job type (`0 = non-production`, `1 = production`)
- `Y`: promotion (`0 = not promoted`, `1 = promoted`)

### Case 2
- `Z`: sex (`0 = female`, `1 = male`)
- `X`: education (`0 = no college`, `1 = college`)
- `Y`: job type (`0 = non-regular`, `1 = regular`)

## Requirements

### R
The R scripts use packages including:
- `dplyr`
- `boot`
- `ggplot2`
- `tikzDevice`
- `lpSolve`
