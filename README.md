# CAMT
Covariate Adaptive Multiple Testing

## Overview
The CAMT package implements a covariate adaptive multiple testing procedure described in the paper: [Covariate Adaptive False Discovery Rate Control with
Applications to Omics-Wide Multiple Testing](https://arxiv.org/***). CAMT allows both the prior null probability and the alternative distribution to depend on covariates. 
It is robust to model mis-specification and is computationally efficient. The package also contains functions for testing the
informativeness of the covariates before using CAMT. It also contains a comprehensive simulation function, which covers all the
settings in the paper.

## Installation         

```
# install.packages("devtools")
devtools::install_github("jchen1981/CAMT")
```



### An Example
We illustrate the usage of CAMT package using simulated data.

```
# Load package
library("CAMT")

# Simulate data
data <- simulate.data(feature.no = 5000, covariate.strength = "Moderate", covariate.model = "pi0",
	sig.density = "Medium", sig.strength = "L3", cor.struct = "None")
  
# Run CAMT  
camt.obj <- camt.fdr(pvals = data$pvals, pi0.var = data$pi0.var, f1.var = data$f1.var, 
	alg.type = "EM", control.method = "knockoff+")
  
# Plot results (logit(pi0) vs covariate, logit(k) vs covariate)
plot.camt(camt.obj, covariate = as.vector(rank(data$pi0.var)), covariate.name = "Covariate rank",
	log = TRUE, file.name = "CovariateModerate")
  
```
