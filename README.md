# CAMT
Covariate Adaptive Multiple Testing v1.1

## Overview
The CAMT package implements two covariate adaptive multiple testing procedures (FDR and FWER) described in [Covariate Adaptive False Discovery Rate Control with
Applications to Omics-Wide Multiple Testing](https://arxiv.org/abs/1909.04811) and [Covariate Adaptive Family-wise Error Control with Applications to Genome-wide Association Studies](https://arxiv.org/abs/2011.01107). CAMT allows  the prior null probability and/or the alternative distribution to depend on covariates. 
It is robust to model mis-specification and is computationally efficient. The package also contains functions for testing the
informativeness of the covariates for multiple testing, and a comprehensive simulation function, which covers a wide range of settings.

## Installation         

```
# install.packages(c("matrixStats", "ggplot2", "cowplot"))
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
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
camt.obj.fdr <- camt.fdr(pvals = data$pvals, pi0.var = data$pi0.var, f1.var = data$f1.var, 
	alg.type = "EM", control.method = "knockoff+")
  
# Plot results (logit(pi0) vs covariate, logit(k) vs covariate)
plot.camt.fdr(camt.obj.fdr, covariate = as.vector(rank(data$pi0.var)), covariate.name = "Covariate rank",
	log = TRUE, file.name = "CovariateModerateFDR.pdf")

camt.obj.fwer <- camt.fwer(pvals = data$pvals, pi0.var = data$pi0.var)
plot.camt.fwer(camt.obj.fwer, covariate = as.vector(rank(data$pi0.var)), covariate.name = 'Covariate rank',
    log = TRUE, file.name = 'CovariateModerateFWER.pdf')	
  
```


