% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/camt.cor.func.R
\name{cov.test.c}
\alias{cov.test.c}
\title{A permutation test for assessing the informativeness of a categorical covariate for CAMT.}
\usage{
cov.test.c(
  pvals,
  covariate,
  cutoffs = quantile(pvals, c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)),
  perm.no = 999,
  n.max = 1e+05,
  silence = TRUE
)
}
\arguments{
\item{pvals}{a numeric vector of p-values.}

\item{covariate}{a factor of the categorical covariate.}

\item{cutoffs}{a numeric vector of the cutoff points for dichotomizing the p-values.}

\item{perm.no}{the number of permutation to assess the significance. Deafult is 999.}

\item{n.max}{an integer number specifying the maximum number of data points to be included. If the number of data points is larger than \code{n.max}, 
subsampling will be performed.}

\item{silence}{a logical value indicating whether to print out the process of the computation.}
}
\value{
A list with the elements
\item{stat.o}{the observed test statistic.}
\item{stat.p}{a vector of the test statistic from permutation.}
\item{p.value}{the omnibu p-value for testing the association.}
\item{p.cut.optim}{the optimal cutoff to dichotomize the p-value.}
}
\description{
The function implements a powerful statistical test for assessing the association between the p-value and a categorical covariate, exploiting the signal sparsity. 
This is achieved by testing the association between two categorical variables after dichotomizing the p-values at the lower end. An omnibus-type
test is designed to combine evidence through various dichotomization. Permutation is used to assess the statistical significance.   We recommend using the 
covariate if the p-value is highly significant (p < 0.005).
}
\examples{

data <- simulate.data(feature.no = 1000, covariate.strength = 'Moderate', covariate.model = 'pi0',
	sig.density = 'Low', sig.strength = 'L3', cor.struct = 'None')
X <- factor(data$pi0.var >= median(data$pi0.var))
obj <- suppressWarnings(cov.test.c(data$pvals, X, perm.no = 2)) 
obj$p.value

}
\references{
Huang JY, ..., Zhang X, Chen J (2020) Leveraging biological and statistical covariates improves the detection power in epigenome-wide
association testing.21(88).
}
\keyword{association}
\keyword{multiple}
\keyword{test}
\keyword{testing,}
