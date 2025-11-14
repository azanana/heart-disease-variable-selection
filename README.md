# Variable Selection Methods for Heart Disease Prediction

Implementation of frequentist and Bayesian variable selection methods for logistic regression, applied to the Cleveland Heart Disease dataset.

## Methods

### Bootstrap + BIC (`bootstrap_bic_variable_selection.R`)
Frequentist approach using bootstrap resampling with BIC-based backward stepwise selection. Computes variable inclusion frequencies across 300 bootstrap samples.

### RJMCMC (`rjmcmc_variable_selection.R`)
Bayesian approach using Reversible Jump Markov Chain Monte Carlo for variable selection. Computes posterior inclusion probabilities for predictors.

## Requirements

```r
library(stats)
library(parallel)
```

## Data

Both scripts expect the Cleveland Heart Disease dataset at `data/processed.cleveland.data`.

Dataset source: UCI Machine Learning Repository  
Citation: Detrano et al. (1989). International application of a new probability algorithm for the diagnosis of coronary artery disease. American Journal of Cardiology, 64(5), 304-310.

## Usage

```r
# Bootstrap + BIC
source("bootstrap_bic_variable_selection.R")

# RJMCMC
source("rjmcmc_variable_selection.R")
```

## Author

Avram Zanana  
Bachelor's Project Mathematics  
University of Groningen, November 2025

## License

MIT License

