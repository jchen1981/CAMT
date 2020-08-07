

process.pvals <- function (pvals) {
	# Process the pvalues to replace 0 or 1
	# Right end 
	pvals[pvals == 1] <- seq(max(pvals[pvals != 1]), 1, len=sum(pvals == 1) + 2)[2 : (sum(pvals == 1) + 1)]
	# Left end
	pvals[pvals == 0]  <- min(pvals[pvals != 0])
	return(pvals)
}

solvek.fdr <- function (pi0, pvals, pvals.cutoff = 1e-15) {
	pvals [pvals  <= pvals.cutoff] <- pvals.cutoff
	pvals [pvals >= 1 - pvals.cutoff] <- 1 - pvals.cutoff
	
	pi1 <- 1 - pi0
	
	fun <- function (k) {
		pk <- pvals ^ (-k)
		sum(- pi1 * pk * ( 1 + (1 - k) * log(pvals)) / (pi0 + pi1 * (1 - k) * pk)) 
	}
	
	err <- try({
				k <- uniroot(fun, c(0.01, 0.99))$root
			}, silent=TRUE)
	
	if (class(err) == 'try-error') {
		k <- 0.75
	} 
	return(k)
}

calculate.weights <- function (pvals, pi0.var = NULL, f1.var = NULL, data = NULL,
		pvals.cutoff = 1e-15, alg.type = c('EM', 'OS'), 
		EM.paras = list(iterlim = 50, tol = 1e-5, k.init = NULL, pi0.init = NULL, nlm.iter = 5),  
		trace = FALSE, return.model.matrix = TRUE, ...) {
	
#		 Compute the covariate-based pi0/k estimates
#		 
#		 Args:
#			 pvals: a numeric vector of p-values.
#			 pi0.var: a formula, a vector, a data frame, or a matrix of covariate values for the prior null probability.
#			 f1.var: a formula, a vector, a data frame, or a matrix of covariate values for the alternative distribution.
#			 data:  a data frame containing the covariates, only used when pi0.var, f1.var are classes of "formula".
#			 pvals.cutoff: a numeric value to replace p-values below that value, which is used to increase the stability of the algorithm.
#			 alg.type:  a character string indicating the algorithm used:
#					 'OS' - direct one-step optimization using a Newton-type algorithm,  
#			         'EM' - using EM algorithm. 
#			         'OS' is much faster but could be inaccurate under some scenarios, e.g., the empirical null is not uniformly distributed.
#			 EM.paras:  a list of control arguments for the EM algorithm:
#			 iterlim - an integer value indicating the maximum number of iterations, 
#			 tol -  a numerical value giving the tolerance in the relative change in the log likelihood below which 
#			 the algorithm is considered to be converged,
#			 pi0.init, k.init - two positive scalars giving the initial guess of the average pi0 and k parameter,
#			 nlm.iter - an integer value indicating the maximal iteration number in running "nlm", used to speed up computation.
#			 trace: a logical value indicating whether to print the process.
#			 return.model.matrix: a logical value indicating whether to return the model matrix. Set FALSE to save memory space.
#			 ...:  parameters to 'nlm' optimization.
#		 
#        Returns:
#		 a list that contains:
#			 call:  the call made.
#			 pi0:   a vector of the estimated null probabilities. 
#			 k:     a vector of the estimated shape parameters for the alternative distribution.
#			 pvals: a numeric vector of the original p-values.
#			 EM.paras:  actually used parameters in EM algorithm.
#			 EM.iter: the iteration used.
#			 loglik:  log likelihood
#			 pi0.coef, k.coef: a vector of the coefficients for pi0 and k.
#			 pi0.var, f1.var: the actual model matrix used if its return is requested.
	alg.type <- match.arg(alg.type)
	
	m0 <- length(pvals)
	nan.ind <- !is.na(pvals)
	pvals <- pvals[nan.ind]
	
	m <- length(pvals)
	
	if (is.null(EM.paras$pi0.init)) {
		if (trace) cat('Estimating initial pi0 values ...\n')
		pi0.est <- qvalue(pvals, pi0.method = 'bootstrap')$pi0
		if (pi0.est >= 0.99)  pi0.est <- 0.99
		EM.paras$pi0.init <- pi0.est
	} else {
		pi0.est <- EM.paras$pi0.init
	}
	
	if (is.null(EM.paras$k.init)) {
		if (trace) cat('Estimating initial k values ...\n')
		k.init <- solvek.fdr(pi0.est, pvals, pvals.cutoff)
		EM.paras$k.init <- k.init
	} else {
		k.init <- EM.paras$k.init
	}
	
	if (is.null(EM.paras$iterlim)) {
		iterlim <- 50
		EM.paras$iterlim <- iterlim
	} else {
		iterlim <- EM.paras$iterlim
	}
	
	if (is.null(EM.paras$nlm.iter)) {
		nlm.iter <- 5
		EM.paras$nlm.iter <- nlm.iter
	} else {
		nlm.iter <- EM.paras$nlm.iter
	}
	
	if (is.null(EM.paras$tol)) {
		tol <- 1e-5
		EM.paras$tol <- tol
	} else {
		tol <- EM.paras$tol
	}
	
	# Threshold the p-value for more robustness
	pvals[pvals <= pvals.cutoff] <- pvals.cutoff
	pvals[pvals >= 1 - pvals.cutoff] <- 1 - pvals.cutoff
	
	if (is.null(pi0.var)) {
		pi0.var <- matrix(rep(1, m))
	} else {
		if (length(intersect(class(pi0.var), c('matrix', 'numeric', 'character', 'factor', 'data.frame', 'formula'))) == 0) {
			stop('Currently do not support the class of pi0.var!\n')
		}
		if (length(intersect('matrix', class(pi0.var))) != 0) {
			if (mode(pi0.var) == 'numeric') {
				if (sum(pi0.var[, 1] == 1) == m0) {
					pi0.var <- pi0.var
				} else {
					pi0.var <- cbind(rep(1, m0), pi0.var)
				}
				
			} else {
				warning('The "pi0.var" matrix contains non-numeric values! Will treat it as a data frame!\n')
				pi0.var <- model.matrix(~., data.frame(pi0.var))
			}
		} else {
			if (length(intersect(c('character', 'factor'), class(pi0.var))) != 0 ) {
				pi0.var <- model.matrix(~., data.frame(pi0.var))
			}  else {
				if (length(intersect(c('data.frame'), class(pi0.var))) != 0) {
					pi0.var <- model.matrix(~., pi0.var)
				}  else {
					if (length(intersect('numeric', class(pi0.var))) != 0) {
						pi0.var <- cbind(rep(1, m0), pi0.var)
						
					} else {
						if (length(intersect(c('formula'), class(pi0.var))) != 0) {
							if (is.null(data)) {
								stop('"pi0.var" is a formula and "data" should be a data frame!\n')
							} else {
								pi0.var <- model.matrix(pi0.var, data)
							}
							
						}
					}
				}
			}
		}
		
		pi0.var <- pi0.var[nan.ind, , drop = FALSE]
	}
	
	
	if (is.null(f1.var)) {
		f1.var <- matrix(rep(1, m))
	} else {
		if (length(intersect(class(f1.var), c('matrix', 'numeric', 'character', 'factor', 'data.frame', 'formula'))) == 0) {
			stop('Currently do not support the class of f1.var!\n')
		}
		if (length(intersect('matrix', class(f1.var))) != 0) {
			if (mode(f1.var) == 'numeric') {
				if (sum(f1.var[, 1] == 1) == m0) {
					f1.var <- f1.var
				} else {
					f1.var <- cbind(rep(1, m0), f1.var)
				}
				
			} else {
				warning('The "f1.var" matrix contains non-numeric values! Will treat it as a data frame!\n')
				f1.var <- model.matrix(~., data.frame(f1.var))
			}
		} else {
			if (length(intersect(c('character', 'factor'), class(f1.var))) != 0 ) {
				f1.var <- model.matrix(~., data.frame(f1.var))
			}  else {
				if (length(intersect(c('data.frame'), class(f1.var))) != 0) {
					f1.var <- model.matrix(~., f1.var)
				}  else {
					if (length(intersect('numeric', class(f1.var))) != 0) {
						f1.var <- cbind(rep(1, m0), f1.var)
						
					} else {
						if (length(intersect(c('formula'), class(f1.var))) != 0) {
							if (is.null(data)) {
								stop('"f1.var" is a formula and "data" should be a data frame!\n')
							} else {
								f1.var <- model.matrix(f1.var, data)
							}
							
						}
					}
				}
			}
		}
		
		f1.var <- f1.var[nan.ind, , drop = FALSE]
	}
	
	
	
	len1 <- ncol(pi0.var)
	len2 <- ncol(f1.var)
	
	if (alg.type %in% c('OS')) {	
		if (trace) cat('Run OS algorithm...\n')
		
		st <- c(binomial()$linkfun(pi0.est), rep(0, ncol(pi0.var) - 1), binomial()$linkfun(k.init), rep(0, ncol(f1.var) - 1))
		
		func1 <- function (x, len1, len2) {
			theta <- x[1 : len1]
			beta <- x[(len1 + 1) : (len1 + len2)]
			exp.eta.theta <- exp(as.vector(pi0.var %*% theta))
			pi0.temp <- (exp.eta.theta / (1 + exp.eta.theta)) 
			exp.eta.beta <- exp(as.vector(f1.var %*% beta))
			k.temp <- (exp.eta.beta / (1 + exp.eta.beta)) 
			res <- -sum(log(pi0.temp  + (1 - pi0.temp) * (1 - k.temp) * pvals  ^ (-k.temp)))
			return(res)	
		}
		suppressWarnings(obj <- nlm(func1, st, len1 = len1, len2 = len2, ...))
		loglik1 <- -obj$minimum
		theta <- obj$estimate[1 : len1]
		beta <- obj$estimate[(len1 + 1) : (len1 + len2)]
		exp.eta.theta <- exp(as.vector(pi0.var %*% theta))
		exp.eta.beta <- exp(as.vector(f1.var %*% beta))
		pi0 <- exp.eta.theta / (1 + exp.eta.theta)
		k <- (exp.eta.beta / (1 + exp.eta.beta))
		
		theta1 <- theta
		beta1 <- beta
	}
	
	if (alg.type %in% c('EM')) {
		
		# No covariate effect - initialization
		theta1 <- c(binomial()$linkfun(pi0.est), rep(0, ncol(pi0.var) - 1))
		beta1 <- c(binomial()$linkfun(k.init), rep(0, ncol(f1.var) - 1))
		
		iter <- 0
		lpval <- log(pvals)
		
		func2 <- function (x, len1, len2, q0) {
			theta <- x[1 : len1]
			beta <- x[(len1 + 1) : (len1 + len2)]
			q1 <- 1 - q0
			exp.eta.theta <- exp(as.vector(pi0.var %*% theta))
			pi0.temp <- (exp.eta.theta / (1 + exp.eta.theta)) 
			exp.eta.beta <- exp(as.vector(f1.var %*% beta))
			k.temp <- (exp.eta.beta / (1 + exp.eta.beta)) 
			res <- -sum(q0 * log(pi0.temp) + q1 * log(1 - pi0.temp) + q1 * (-k.temp * lpval + log(1 - k.temp))) 
			return(res)
		}
		
		loglik0 <- 0
		if (trace) cat('Run EM algorithm...\n')
		while (iter < iterlim) {
			if (trace) cat('.')
			# E step
			exp.eta.theta <- exp(as.vector(pi0.var %*% theta1))
			exp.eta.beta <- exp(as.vector(f1.var %*% beta1))
			pi0 <- exp.eta.theta / (1 + exp.eta.theta)
			k <- exp.eta.beta / (1 + exp.eta.beta)
			f1 <- (1 - k) * pvals  ^ (-k)
			f01 <-  (1 - pi0) * f1 + pi0
			q0 <- pi0  / f01
			
			# M step
			st <- c(theta1, beta1)
			obj <- suppressWarnings(nlm(func2, st, len1 = len1, len2 = len2, q0 = q0,  iterlim = nlm.iter, ...))
			
			# Update new likelihood
			theta2 <- obj$estimate[1 : len1]
			beta2 <- obj$estimate[(len1 + 1) : (len1 + len2)]
			
			loglik1 <- sum(log(f01))
			
			if (abs((loglik1 - loglik0) / loglik0) < tol) {
				break
			} else {
				theta1 <- theta2
				beta1 <- beta2
				iter <- iter + 1
				loglik0 <- loglik1
			}
		}
		
#		if (iter == iterlim) warning('Maximum number of iterations reached!')
		if (trace) cat('\n')
	}
	
	temp <- rep(NA, m)
	temp[nan.ind] <- pi0
	pi0 <- temp
	
	result <- list(call = match.call(),  pi0 = pi0,  k = k, pi0.coef = theta1, k.coef = beta1, 
			EM.paras = EM.paras,  loglik = loglik1, pvals = pvals)
	
	if (return.model.matrix) {
		result$pi0.var <- result$pi0.var
		result$f1.var <- result$f1.var
	}
	
	if (alg.type == 'EM') {
		result$EM.iter <- iter
	}
	return(result)
	
}



hybridfdr.full.path <- function (pvals, pi0, k,  control.method=c('hybrid', 'knockoff+'), burnin.no = 500,
		pi0.low = 0.1, pi0.high = 0.9999) {
	
#		 Perform FDR control based on the estimated pi0/k.
#		 
#		 Args:
#		     pvals: a numeric vector of p-values.
#			 pi0:   a numeric vector of the estimated null probabilities. 
#			 k:     a numeric vector of the estimated shape parameters for the alternative distribution.
#			 control.method: a character string indicating the FDR control variant:
#					 knockoff+ - the knockoff+ procedure of Barber-Candes (BC), conservative at sparse signal/small FDR levels
#			         hybrid -  an empirical method to correct the conservativeness of 'knockoff+'. The method is
#			            based on taking the maximum over the BC-type (knockoff) and BH-type FDR estimates (hybrid method) for a certain number
#			            ("burnin.no") of the most promising hypotheses at the start of the algorithm. The rest use knockoff-based FDR estimator.
#			 burnin.no: an integer value indicating the number of the most promising hypotheses that will apply the hybrid procedure above.
#			 pi0.low, pi0.high: the lowest/highest allowed pi0 values, which could guard against the dominance of the prior.
#		 
#		 Returns:
#		 a list that contains:
#			 fdr:   a numeric vector of the adjusted p-values. 
#			 ts:    a numeric vector of the thresholds (t) below which the corresponding hypothesis will be rejected.
#			 pi0:   a numeric vector of the actual pi0 used (truncated at pi0.low)
	
	control.method <- match.arg(control.method)
	
	pvals <- process.pvals(pvals)
	
	pi0[pi0 < pi0.low] <- pi0.low
	pi0[pi0 == 1] <- pi0.high
	
	n <- length(pvals)
	rs <-  (1 - pi0) / pi0
	#rs <- 1 / rs
	qs <- dbeta(pvals, shape1 = 1 - k, shape2 = 1) * rs
	ts <- 1 / (1 + qs)
	
	out <- sort(ts,  index.return = TRUE)
	ts <- out$x
	index0 <- out$ix
	
	rs <- rs[index0]
	qs <- qs[index0]
	pvals <- pvals[index0]
	k <- k[index0]
	pi0 <- pi0[index0]
	
	# ts cutoffs for all the p-values
	ts.p <- 1 / (1 + (1 - k) * (1 - pvals) ^ (-k) * rs)
	
	ts.p <- sort(c(ts, ts.p), decreasing = TRUE)
	m <- length(ts.p)
	
	if (control.method == 'knockoff+') {
		fdr <-  (1 + (m + 1) - match(ts, ts.p) - 1:n) / (1:n)
	} 
	
	if (control.method == 'hybrid') {
		fdr <-  (0 + (m + 1) - match(ts, ts.p) - 1:n) / (1:n)
		burnin.no <- min(burnin.no, length(pi0))
		fdr.burnin <- sapply(1:burnin.no, function (i) {
					t0 <- ts[i]
					cs <- pmax((1 - t0) / t0 / rs, 1 - k)
					sum(pi0 * (cs / (1 - k)) ^ (-1 / k)) / i 
				})
		fdr[1:burnin.no] <- pmax(fdr[1:burnin.no], fdr.burnin)
	}
	
	fdr <- rev(pmin(1, cummin(rev(fdr))))
	index0 <- order(index0)
	
	return(list(fdr = fdr[index0], ts = ts[index0], pi0 = pi0[index0]))
}


#' Perform covariate-adaptive false discovery rate control
#'
#' The function implements a scalable, flexible, robust and powerful FDR control method for large-scale multiple testing exploiting the auxiliary covariates. 
#' It allows both the prior null probability and the alternative distribution to depend on covariates.
#'
#' @param pvals a numeric vector of p-values.
#' @param pi0.var a formula, a vector, a data frame, or a matrix of covariate values for the prior null probability.
#' @param f1.var a formula, a vector, a data frame, or a matrix of covariate values for the alternative distribution.
#' @param data  a data frame containing the covariates, only used when pi0.var, f1.var are classes of 'formula'.
#' @param pvals.cutoff a numeric value to replace p-values below that value, which is used to increase the stability of the algorithm. 
#' @param pi0.low the allowed minimal pi0 value, which could guard against the dominance of the prior. 
#' @param alg.type  a character string indicating the algorithm used. 'OS' - direct one-step optimization using a Newton-type algorithm, 
#' 'EM' - using EM algorithm. 'OS' is fast but could be inaccurate under some scenarios. Default 'EM'.
#' @param EM.paras  a list of control arguments for the EM algorithm
#' \itemize{
#' \item{iterlim}{an integer value indicating the maximum number of iterations.}
#' \item{tol}{a numeric value giving the tolerance in the relative change in the log likelihood below which the algorithm is considered to be converged.}
#' \item{pi0.init, k.init}{two scalars giving the initial guess of the average pi0 and k parameter.}
#' \item{nlm.iter}{an integer indicating the allowed maximum iteration in running \code{'nlm'}, used to speed up computation.}
#' }
#' @param control.method  a character string indicating the FDR control variant. 'knockoff+': the knockoff+ procedure of Barber-Candes (BC), 
#' conservative at sparse signal/small FDR levels. 'hybrid': an empirical method to correct the conservativeness of 'knockoff+'. The method is 
#' based on taking the maximum over the BC-type (knockoff) and BH-type FDR estimates for a certain number (as specified by \code{'burnin.no'}) of 
#' the most promising hypotheses at the start of the algorithm. The rest use knockoff-based FDR estimator. Default is 'hybrid'.
#' @param burnin.no  an integer value indicating the number of the most promising hypotheses that will apply the 'hybrid' procedure above.
#' @param trace a logical value indicating whether to print the process.
#' @param return.model.matrix a logical value indicating whether to return the model matrix.  Consider setting to FALSE if it's huge.
#' @param ... parameters passing to \code{'nlm'} optimization.
#'
#' @return A list with the elements
#' \item{call}{the call made.}
#' \item{pi0}{a vector of the estimated null probabilities.}
#' \item{k}{a vector of the estimated shape parameters for the alternative distribution.}
#' \item{EM.paras}{actually used parameters in EM algorithm.}
#' \item{EM.iter}{the number of iteration actually used.}
#' \item{loglik}{log likelihood.}
#' \item{pi0.coef, k.coef}{a vector of the coefficients for pi0 and k.}
#' \item{pi0.var, f1.var}{the actual model matrix used if its return is requested.}
#' \item{fdr}{a numeric vector of the adjusted p-values.}
#' \item{pvals}{a numeric vector of the original p-values used.}
#' \item{ts}{a numeric vector of the thresholds (t) below which the corresponding hypothesis will be rejected.}
#' 
#' @author Jun Chen
#' @references Xianyang Zhang, Jun Chen. Covariate Adaptive False Discovery Rate Control with Applications to Omics-Wide Multiple Testing. JASA. To appear.
#' @keywords FDR
#' @importFrom stats binomial dbeta qnorm uniroot nlm model.matrix
#' @importFrom qvalue qvalue
#' @examples
#'
#' data <- simulate.data(feature.no = 10000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Medium', sig.strength = 'L3', cor.struct = 'None')
#' camt.fdr.obj <- camt.fdr(pvals = data$pvals, pi0.var = data$pi0.var, f1.var = data$f1.var, 
#'		alg.type = 'EM', control.method = 'knockoff+')
#' plot.camt.fdr(camt.fdr.obj, covariate = as.vector(rank(data$pi0.var)), covariate.name = 'Covariate rank',
#'		log = TRUE, file.name = 'CovariateModerate.pdf')
#' 
#' @rdname camt.fdr
#' @export

camt.fdr <- function (pvals, pi0.var = NULL, f1.var = NULL, data = NULL,
		pvals.cutoff = 1e-15, pi0.low = 0.1,  alg.type = c('EM', 'OS'), 
		EM.paras = list(iterlim = 50, k.init = NULL, tol = 1e-5, pi0.init = NULL, nlm.iter = 5),  
		control.method = c('hybrid', 'knockoff+'), burnin.no = 500,
		trace = FALSE, return.model.matrix = TRUE, ...) {
	

	control.method <- match.arg(control.method)
	alg.type <- match.arg(alg.type)
	
	if (length(pvals) <= 2000) {
		warning(paste0("The number of feature size is less than 2,000. The procedure will be VERY conservative unless the covariate is highly informative ",
						"or the signal is dense! Please first use 'cov.test.n' or 'cov.test.c'to test the informativeness of the covariate.\n"))
	} 
	
	if (trace) cat('Calculate genomic inflation factor ...\n')
	lambda <- median(qchisq(pvals[pvals >= 0.5], df = 1, lower.tail = FALSE)) / 0.102
	
	if (trace) cat('Genomic inflation factor for p-value [0.5, 1] is ', lambda, '\n')
	
	if (lambda >= 1.1) {
		warning("The p-value distribution under the null seems right-skewed! The procedure will be anticonservative!\n")
	} else {
		if (lambda <= 0.9) {
			warning("The p-value distribution under the null seems left-skewed! The procedure will be conservative!\n")
		} else {
			if (trace) cat("The p-value distribution under the null seems OK!\n")
		}
	}
	
	
	if (trace) cat('Estimate pi0 and k ...\n')
	camt.fdr.obj <- calculate.weights(pvals = pvals, pi0.var = pi0.var, f1.var = f1.var, data = data,
			pvals.cutoff = pvals.cutoff, alg.type = alg.type, 
			EM.paras = EM.paras, trace = trace, ...)
	
	if (trace) cat('Perform covariate-adaptive FDR control ...\n')
	fdr.obj <- hybridfdr.full.path(pvals, camt.fdr.obj$pi0, camt.fdr.obj$k, pi0.low = pi0.low, control.method = control.method, burnin.no = burnin.no)
	
	if (trace) cat('Finished!\n')
	result <- camt.fdr.obj
	result$fdr <- fdr.obj$fdr
	result$ts <- fdr.obj$ts
	result$pi0 <- fdr.obj$pi0
	
	return(result)
	
}



#' Plot the camt.fdr results
#'
#' The function plots the null probabilities, f1 shape parameters and p-values against covariate values. Significant hypotheses are highlighted in red.
#'
#' @param camt.fdr.obj a list returned from running 'camt.fdr'.
#' @param covariate a vector containing the one-dimensional covariate values.
#' @param covariate.type a character string of either "continuous" or "categorical".
#' @param covariate.name a character string for the name of the covariate to be plotted.
#' @param alpha  a numeric value, the target FDR level.
#' @param nlimit an integer, the number of insignificant hypotheses to be sampled to reduce the visualization complexity.
#' @param log a logical value indicating whether the p-values should be plotted on the log scale. 
#' @param logit a logical value indicating whether the pi0/k should be plotted on the logit scale. 
#' @param file.name file.name a character string (ended with .pdf) for the name of the generated pdf file. Could include the path. 
#'
#' @return A list of \code{'ggplot2'} objects.
#' 
#' @author Jun Chen
#' @references Xianyang Zhang, Jun Chen. Covariate Adaptive False Discovery Rate Control with Applications to Omics-Wide Multiple Testing. JASA. To appear.
#' @keywords FDR, multiple testing
#' @import ggplot2
#' @import grDevices
#' @importFrom cowplot plot_grid
#' @examples
#'
#' data <- simulate.data(feature.no = 10000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Medium', sig.strength = 'L3', cor.struct = 'None')
#' camt.fdr.obj <- camt.fdr(pvals = data$pvals, pi0.var = data$pi0.var, f1.var = data$f1.var, 
#'		alg.type = 'EM', control.method = 'knockoff+')
#' plot.camt.fdr(camt.fdr.obj, covariate = as.vector(rank(data$pi0.var)), covariate.name = 'Covariate rank',
#'		log = TRUE, file.name = 'CovariateModerate.pdf')
#' 
#' @rdname plot.camt.fdr
#' @export

plot.camt.fdr <- function (camt.fdr.obj,  covariate, covariate.type = c('continuous',  'categorical'), covariate.name = 'covariate', 
		alpha = 0.1, nlimit = 10000, log = FALSE, logit = TRUE,  file.name = 'CAMT.plot.pdf') {
	
	covariate.type <- match.arg(covariate.type)
	pvals <- camt.fdr.obj$pvals
	pi0 <- camt.fdr.obj$pi0
	k <- camt.fdr.obj$k
	fdr <- camt.fdr.obj$fdr
	ts <- camt.fdr.obj$ts
	
	# remove NA
	ind <- !is.na(pvals) 
	pvals <- pvals[ind]
	pi0 <- pi0[ind]
	k <- k[ind]
	fdr <- fdr[ind]
	covariate <- covariate[ind]
	ts <- ts[ind]
	
	n <- length(pvals)
	
	rejection <- fdr <= alpha
	cat(paste0('CAMT identified ', sum(rejection), ' at an FDR of ', alpha, '.\n'))
	
	ind1 <- which(rejection)
	ind2 <- setdiff(1:n, ind1)
	
	if (length(ind2) > nlimit) {
		ind2 <- sample(ind2, nlimit)
	}
	
	ind <- c(ind1, ind2)
	rejection <- factor(c(rep(1, length(ind1)), rep(0, length(ind2))))
	
	pvals <- pvals[ind]
	pi0 <- pi0[ind]
	k <- k[ind]
	fdr <- fdr[ind]
	covariate <- covariate[ind]
	ts <- ts[ind]
	
	lfdr <- (pi0 /  (pi0 + (1 - pi0) * dbeta(pvals, shape1 = 1 - k, shape2 = 1)))
	
	pi0[pi0 == 0] <- min(pi0[pi0 != 0])
	pi0[pi0 == 1] <- max(pi0[pi0 != 1])
	
	if (log == TRUE) {
		pvals <- -log10(pvals)
		ylab.pval <- '-log10(p-value)'
	} else {
		ylab.pval <- 'p-value'
	}
	
	if (logit == TRUE) {
		pi0 <- log(pi0 / (1 - pi0))
		k <- log(k / (1 - k))
		ylab.pi0 <- 'logit(pi0)'
		ylab.k <- 'logit(k)'
	} else {
		ylab.pi0 <- 'pi0'
		ylab.k <- 'k'
	}
	
	data <- data.frame(pi0 = pi0, k = k, fdr = fdr, lfdr = lfdr, rejection = rejection, covariate = covariate, pvals = pvals)
	data <- data[rev(order(data$fdr)), ]
	
	plot.list <- list()
	if (covariate.type == 'continuous') {
		data1 <- subset(data, rejection == '0')
		data2 <- subset(data, rejection == '1')
		plot.list[['pi0']] <- ggplot(data, aes(x = covariate, y = pi0, col = rejection)) + 
				geom_point(aes(x = covariate, y = pi0), data = data1, col = 'darkgray', alpha = 0.75, size = 0.25) +
				geom_point(aes(x = covariate, y = pi0), data = data2, col = 'red', alpha = 0.75, size = 0.25) +
				ylab(ylab.pi0) +
				xlab(covariate.name) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				ggtitle('Null probability vs covariate') +
				theme_bw() +
				theme(legend.position = "none")
		
		plot.list[['k']] <- ggplot(data, aes(x = covariate, y = k, col = rejection)) + 
				geom_point(alpha = 0.75, size = 0.25) +
				#	geom_smooth(aes(x = covariate, y = k), method = 'lm', inherit.aes = FALSE, size = 0.25) +
				ylab(ylab.k) +
				ggtitle('f1 shape parameter vs covariate') +
				scale_colour_manual(values = c('darkgray', 'red')) +
				xlab(covariate.name) +
				theme_bw()
		
		
		plot.list[['pval']] <- ggplot(data, aes(x = covariate, y = pvals, col = rejection)) + 
				geom_point(alpha = 0.75, size = 0.2) +
				#	geom_area(aes(x = covariate, y = pvals.c), col = 'red', alpha = 0.25, size = 0.25) +
				ylab(ylab.pval) +
				ggtitle(paste0('Target FDR level (alpha = ', alpha, ')')) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				xlab(covariate.name) +
				theme_bw() +
				theme(legend.position = "none")
		
		plot.list[['lfdr']] <- ggplot(data, aes(x = covariate, y = pvals, col = lfdr)) + 
				geom_point(alpha = 0.75, size = 0.2) +
				ylab(ylab.pval) +
				ggtitle('Estimated local FDR') +
				xlab(covariate.name) +
				scale_color_gradient(low="red", high="darkgray") +
				theme_bw()
	}
	
	if (covariate.type == 'categorical') {
		plot.list[['pi0']] <- ggplot(data, aes(x = covariate, y = pi0, col = rejection)) + 
				geom_boxplot(alpha = 0.75) +
				#	geom_smooth(aes(x = covariate, y = pi0), method = 'lm', inherit.aes = FALSE, size = 0.25) +
				ylab(ylab.pi0) +
				xlab(covariate.name) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				ggtitle('Null probability vs covariate') +
				theme_bw() +
				theme(legend.position = "none")
		
		plot.list[['k']] <- ggplot(data, aes(x = covariate, y = k, col = rejection)) + 
				geom_boxplot(alpha = 0.75) +
				#	geom_smooth(aes(x = covariate, y = k), method = 'lm', inherit.aes = FALSE, size = 0.25) +
				ylab(ylab.k) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				ggtitle('f1 shape parameter vs covariate') +
				xlab(covariate.name) +
				theme_bw()
		
		
		plot.list[['pval']] <- ggplot(data, aes(x = covariate, y = pvals, col = rejection)) + 
				#	geom_boxplot() +
				geom_jitter(alpha = 0.75, size = 0.2, position=position_jitter(height=0)) +
				#	geom_area(aes(x = covariate, y = pvals.c), col = 'red', alpha = 0.25, size = 0.25) +
				ylab(ylab.pval) +
				ggtitle(paste0('Target FDR level (alpha = ', alpha, ')')) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				xlab(covariate.name) +
				theme_bw() +
				theme(legend.position = "none")
		
		plot.list[['lfdr']] <- ggplot(data, aes(x = covariate, y = pvals, col = lfdr)) + 
				#	geom_boxplot() +
				geom_jitter(alpha = 0.75, size = 0.2, position=position_jitter(height=0)) +
				ylab(ylab.pval) +
				ggtitle('Estimated local FDR') +
				xlab(covariate.name) +
				scale_color_gradient(low="red", high="darkgray") +
				theme_bw()
	}
	cat('Generate combined PDF plot ...\n')
	pdf(paste0(file.name, '.pdf'), width =10, height = 6)
	obj <- plot_grid(plot.list[['pi0']], plot.list[['k']], plot.list[['pval']], plot.list[['lfdr']], ncol=2)
	print(obj)
	dev.off()
	
	cat('Finished!\n')
	return(invisible(plot.list))
}

############################################################################################################
# FWER procedure
solvek.fwer <- function(pvals, pi0, init = 0.25, nlm.iter = 10, pvals.cutoff = 1e-15) {
	pvals [pvals  <= pvals.cutoff] <- pvals.cutoff
	pvals [pvals >= 1 - pvals.cutoff] <- 1 - pvals.cutoff
	fun <- function(k) {
		- sum(log(pi0 + (1 - pi0) * k * pvals ^ (k - 1))) 
	}
	
	err <- try(k <- nlminb(init, fun, control = list(iter.max = nlm.iter), 
					lower = 0, upper = 1)$par, silent = TRUE)
	if (class(err) == 'try-error') {
		k <- 0.25
	} 
	return(k)
}

#' Return the rejected hypotheses at a given significance level.
#'
#' The function outputs the indices of the rejected hypotheses at a specified signficance level.
#'
#' @param camt.fwer.obj the output from running \code{'camt.fwer'}.
#' @param alpha a numeric value for the desired signficance level.
#'
#' @return A vector containing the indices of the rejected hypotheses.
#' 
#' @author Huijuan Zhou
#' @references Huijuan Zhou, Xianyang Zhang, Jun Chen. Covariate Adaptive Family-wise Error Control with Applications to Genome-wide Association Studies. Submitted.
#' @keywords FWER, multiple testing
#' @examples
#'
#' data <- simulate.data(feature.no = 10000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Medium', sig.strength = 'L4', cor.struct = 'None')
#' camt.obj.fwer <- camt.fwer(pvals = data$pvals, pi0.var = data$pi0.var)
#' camt.fwer.rejection(camt.obj.fwer, alpha = 0.1)
#' 
#' @rdname camt.fwer.rejection
#' @export
camt.fwer.rejection <- function (camt.fwer.obj, alpha = 0.1) {
	
	pvals <- camt.fwer.obj$pvals
	pi0 <- camt.fwer.obj$pi0
	gamma <- camt.fwer.obj$gamma
	k <- camt.fwer.obj$k
	Y <- (pvals > gamma)
	
	pi1 <- 1 - pi0
	tau <- k * (alpha * (1 - gamma)) ^ (k - 1) * (sum(Y * (pi1 / pi0) ^ (1 / (1 - k)))) ^ (1 - k)
	t <- (pi1 * k / (pi0 * tau)) ^ (1 / (1 - k))
	rej <- which(pvals < pmin(t, gamma))
	return(rej)
}

#' Perform the covariate-adaptive family-wise error rate control
#'
#' The function implements a scalable, flexible, robust and powerful FWER control method for large-scale multiple testing exploiting the auxiliary covariates. 
#' It allows both the prior null probability to depend on covariates.
#'
#' @param pvals a numeric vector of p-values.
#' @param pi0.var a formula, a vector, a data frame, or a matrix of covariate values for the prior null probability.
#' @param data  a data frame containing the covariates, only used when pi0.var are classes of 'formula'.
#' @param pvals.cutoff a numeric value to replace p-values below that value, which is used to increase the stability of the algorithm. 
#' @param alpha a numeric value for the significance level. Default is 0.1.
#' @param k.fix a numeric value between 0 and 1. It is the k parameter of the alternative (beta) distribution.
#' If it is null, it will be estimated via EM.
#' @param EM.paras  a list of control arguments for the EM algorithm
#' \itemize{
#' \item{iterlim}{an integer value indicating the maximum number of iterations.}
#' \item{tol}{a numeric value giving the tolerance in the relative change in the log likelihood below which the algorithm is considered to be converged.}
#' \item{pi0.init, k.init}{two scalars giving the initial guess of the  pi0 and k parameter.}
#' \item{nlm.iter}{an integer indicating the allowed maximum iteration in running \code{'nlm'}, used to speed up computation.}
#' }
#' @param trace a logical value indicating whether to print the process.
#' @param return.model.matrix a logical value indicating whether to return the model matrix.  Consider setting to FALSE if it's huge.
#'
#' @return A list with the elements
#' \item{call}{the call made.}
#' \item{pi0}{a vector of the estimated null probabilities.}
#' \item{k}{the (estimated) value for the k parameter of the alternative distribution.}
#' \item{pi0.coef}{a vector of the coefficients for pi0.}
#' \item{gamma}{the estimated gamma value to dichotomize the p-values}
#' \item{loglik}{log likelihood.}
#' \item{EM.paras}{actually used parameters in EM algorithm.}
#' \item{EM.iter}{the number of iteration actually used.}
#' \item{rejection}{a vector containing the indices of the rejected hypotheses.}
#' 
#' @author Huijuan Zhou
#' @references Huijuan Zhou, Xianyang Zhang, Jun Chen. Covariate Adaptive Family-wise Error Control with Applications to Genome-wide Association Studies. Submitted.
#' @keywords FWER
#' @importFrom stats binomial dbeta qnorm uniroot nlm nlminb
#' @importFrom qvalue qvalue
#' @examples
#'
#' data <- simulate.data(feature.no = 10000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Medium', sig.strength = 'L4', cor.struct = 'None')
#' camt.obj.fwer <- camt.fwer(pvals = data$pvals, pi0.var = data$pi0.var)
#' plot.camt.fwer(camt.obj.fwer, covariate = as.vector(rank(data$pi0.var)), covariate.name = 'Covariate rank',
#'		log = TRUE, file.name = 'CovariateModerate.pdf')
#' 
#' @rdname camt.fwer
#' @export

camt.fwer <- function(pvals, pi0.var, data = data, pvals.cutoff = 1e-15, alpha = 0.1, k.fix = NULL,
		EM.paras = list(iterlim = 50, tol = 1e-5, k.init = NULL, 
				pi0.init = NULL, nlm.iter = 5), trace = FALSE, return.model.matrix = FALSE) {
	
	m0 <- length(pvals)
	nan.ind <- !is.na(pvals)
	pvals <- pvals[nan.ind]
	m <- length(pvals)
	
	pvals[pvals == 1] <- seq(max(pvals[pvals != 1]), 1, 
			length.out = sum(pvals == 1) + 2)[2 : (sum(pvals == 1) + 1)]
	pvals[pvals == 0] <- min(pvals[pvals != 0])
	qvals <- qvalue(pvals, pi0.method = "bootstrap")
	gamma <- qvals$lambda[which.min(abs(qvals$pi0.lambda - qvals$pi0))]
	Y <- (pvals > gamma)
	b0 <- (1 - gamma) * Y + gamma * (1 - Y)
	
	# Initialize pi0
	if (is.null(EM.paras$pi0.init)) {
		if (trace) cat('Estimating initial pi0 values ...\n')
		pi0.init <- qvals$pi0
		if (pi0.init >= 0.99)  pi0.init <- 0.99
		EM.paras$pi0.init <- pi0.init
	} else {
		pi0.init <- EM.paras$pi0.init
	}
	
	# Initialize k
	if (is.null(EM.paras$k.init)) {
		if (trace) cat('Estimating initial k values ...\n')
		k.init <- solvek.fwer(pvals, pi0.init, pvals.cutoff = pvals.cutoff)
		EM.paras$k.init <- k.init
	} else {
		k.init <- EM.paras$k.init
	}
	
	if (is.null(EM.paras$iterlim)) {
		iterlim <- 50
		EM.paras$iterlim <- iterlim
	} else {
		iterlim <- EM.paras$iterlim
	}
	
	if (is.null(EM.paras$nlm.iter)) {
		nlm.iter <- 5
		EM.paras$nlm.iter <- nlm.iter
	} else {
		nlm.iter <- EM.paras$nlm.iter
	}
	
	if (is.null(EM.paras$tol)) {
		tol <- 1e-5
		EM.paras$tol <- tol
	} else {
		tol <- EM.paras$tol
	}
	
	if (is.null(pi0.var)) {
		pi0.var <- matrix(rep(1, m))
	} else {
		if (length(intersect(class(pi0.var), c('matrix', 'numeric', 'character', 'factor', 'data.frame', 'formula'))) == 0) {
			stop('Currently do not support the class of pi0.var!\n')
		}
		if (length(intersect('matrix', class(pi0.var))) != 0) {
			if (mode(pi0.var) == 'numeric') {
				if (sum(pi0.var[, 1] == 1) == m0) {
					pi0.var <- pi0.var
				} else {
					pi0.var <- cbind(rep(1, m0), pi0.var)
				}
				
			} else {
				warning('The "pi0.var" matrix contains non-numeric values! Will treat it as a data frame!\n')
				pi0.var <- model.matrix(~., data.frame(pi0.var))
			}
		} else {
			if (length(intersect(c('character', 'factor'), class(pi0.var))) != 0 ) {
				pi0.var <- model.matrix(~., data.frame(pi0.var))
			}  else {
				if (length(intersect(c('data.frame'), class(pi0.var))) != 0) {
					pi0.var <- model.matrix(~., pi0.var)
				}  else {
					if (length(intersect('numeric', class(pi0.var))) != 0) {
						pi0.var <- cbind(rep(1, m0), pi0.var)
						
					} else {
						if (length(intersect(c('formula'), class(pi0.var))) != 0) {
							if (is.null(data)) {
								stop('"pi0.var" is a formula and "data" should be a data frame!\n')
							} else {
								pi0.var <- model.matrix(pi0.var, data)
							}
							
						}
					}
				}
			}
		}
		
		pi0.var <- pi0.var[nan.ind, , drop = FALSE]
	}
	
	# Threshold the p-value for more robustness
	pvals[pvals <= pvals.cutoff] <- pvals.cutoff
	pvals[pvals >= 1 - pvals.cutoff] <- 1 - pvals.cutoff
	
	len.theta <- ncol(pi0.var)
	theta <- c(binomial()$linkfun(pi0.init), rep(0, len.theta - 1))
	
	if (is.null(k.fix)) {
		k <- k.init
	} else {
		k <- k.fix
	}
	
	fun1 <- function(theta, q0) {
		q1 <- 1 - q0
		exp.eta.theta <- exp(as.vector(pi0.var %*% theta))
		pi0 <- exp.eta.theta / (1 + exp.eta.theta)
		- sum(q0 * log(pi0) + q1 * log(1 - pi0))
	}
	fun2 <- function(k, q0) {
		q1 <- 1 - q0  
		gammak <- gamma ^ k
		- sum(q1 * (Y * log(1 - gammak) + k * (1 - Y) * log(gamma)))
	}
	
	iter <- 1
	loglik1 <- 1
	
	if (trace) cat('Run EM algorithm...\n')
	while (iter <= iterlim) {
		if (trace) cat('.')
		# E step
		exp.eta.theta <- exp(as.vector(pi0.var %*% theta))
		pi0 <- exp.eta.theta / (1 + exp.eta.theta)
		
		gammak <- gamma ^ k
		b1 <- (1 - gammak) * Y + gammak * (1 - Y)
		pi0b0 <- pi0 * b0
		pib <- pi0b0 + (1 - pi0) * b1
		q0 <- pi0b0 / pib
		
		loglik2 <- sum(log(pib))
		
		if ((abs((loglik2 - loglik1) / loglik1) < tol) | (iter >= iterlim)) {
			break
		} else {
			theta <- suppressWarnings(nlm(fun1, theta, q0 = q0, iterlim = nlm.iter))$estimate
			if (is.null(k.fix)) {
				k <- suppressWarnings(nlminb(k, fun2, q0 = q0, control = list(iter.max = nlm.iter), 
								lower = 0, upper = 1))$par
			}
		}
		loglik1 <- loglik2
		iter <- iter + 1
	}
	
	if (iter == iterlim) warning('Maximum number of iterations reached!')
	
	temp <- rep(NA, m)
	temp[nan.ind] <- pi0
	pi0 <- temp
	
	camt.fwer.obj <- list(call = match.call(), pvals = pvals, pi0.var = pi0.var, alpha = alpha,
			pi0 = pi0, k = k, pi0.coef = theta, gamma = gamma,
			loglik = loglik2, EM.paras = EM.paras, EM.iter = iter) 
	
	camt.fwer.obj$rejection <- camt.fwer.rejection(camt.fwer.obj, alpha)
	
	if (return.model.matrix) {
		camt.fwer.obj$pi0.var <- camt.fwer.obj$pi0.var
	}
	
	return(camt.fwer.obj)
}


#' Plot the camt.fwer results
#'
#' The function plots the null probabilities and p-values against covariate values. Significant hypotheses are highlighted in red.
#'
#' @param camt.fwer.obj a list returned from running 'camt.fwer'.
#' @param covariate a vector containing the one-dimensional covariate values.
#' @param covariate.type a character string of either "continuous" or "categorical". 
#' @param covariate.name a character string for the name of the covariate to be plotted.
#' @param alpha  a numeric value, the target FWER level.
#' @param nlimit an integer, the number of insignificant hypotheses to be sampled to reduce the visualization complexity.
#' @param log a logical value indicating whether the p-values should be plotted on the log scale. 
#' @param logit a logical value indicating whether the pi0 should be plotted on the logit scale. 
#' @param file.name a character string (ended with .pdf) for the name of the generated pdf file. Could include the path. 
#'
#' @return A list of \code{'ggplot2'} objects.
#' 
#' @author Huijuan Zhou
#' @references Huijuan Zhou, Xianyang Zhang, Jun Chen. Covariate Adaptive Family-wise Error Control with Applications to Genome-wide Association Studies. Submitted.
#' @keywords FWER, multiple testing
#' @import ggplot2
#' @import grDevices
#' @importFrom cowplot plot_grid
#' @examples
#'
#' data <- simulate.data(feature.no = 10000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Medium', sig.strength = 'L4', cor.struct = 'None')
#' camt.obj.fwer <- camt.fwer(pvals = data$pvals, pi0.var = data$pi0.var)
#' plot.camt.fwer(camt.obj.fwer, covariate = as.vector(rank(data$pi0.var)), covariate.name = 'Covariate rank',
#'		log = TRUE, file.name = 'CovariateModerate.pdf')
#' 
#' @rdname plot.camt.fwer
#' @export
plot.camt.fwer <- function(camt.fwer.obj, covariate, covariate.type = c('continuous',  'categorical'), covariate.name = 'covariate',
		alpha = 0.1, nlimit = 10000, log = FALSE, logit = TRUE, file.name = 'CAMT.fwer.plot') {
	
	covariate.type <- match.arg(covariate.type)
	pvals <- camt.fwer.obj$pvals
	pi0.var <- covariate
#  alpha <- camt.fwer.obj$alpha
	
	pi0 <- camt.fwer.obj$pi0
	
#  rejection <- camt.fwer.obj$rejection
	rejection <- camt.fwer.rejection(camt.fwer.obj, alpha)
	
	cat(paste0('CAMT identified ', length(rejection), ' at an FWER of ', alpha, '.\n'))
	
	m <- length(pvals)
	ind1 <- rejection
	ind2 <- setdiff(1 : m, ind1)
	
	if (length(ind2) > nlimit) {
		ind2 <- sample(ind2, nlimit)
	}
	
	ind <- c(ind1, ind2)
	pvals <- pvals[ind]
	pi0.var <- pi0.var[ind]
	pi0 <- pi0[ind]
	rejection <- factor(c(rep(1, length(ind1)), rep(0, length(ind2))))
	
	if (log == TRUE) {
		pvals <- -log10(pvals)
		ylab.pval <- '-log10(p-value)'
	} else {
		ylab.pval <- 'p-value'
	}
	
	if (logit == TRUE) {
		pi0 <- log(pi0 / (1 - pi0))
		ylab.pi0 <- 'logit(pi0)'
	} else {
		ylab.pi0 <- 'pi0'
	}
	
	data <- data.frame(pi0 = pi0, rejection = rejection, pvals = pvals, pi0.var = pi0.var)
	
	plot.list <- list()
	if (covariate.type == 'continuous') {
		data1 <- subset(data, rejection == '0')
		data2 <- subset(data, rejection == '1')
		plot.list[['pi0']] <- ggplot(data, aes(x = pi0.var, y = pi0, col = rejection)) + 
				geom_point(aes(x = pi0.var, y = pi0), data = data1, col = 'darkgray', alpha = 0.75, size = 0.25) +
				geom_point(aes(x = pi0.var, y = pi0), data = data2, col = 'red', alpha = 0.75, size = 0.25) +
				xlab(covariate.name) +
				ylab(ylab.pi0) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				ggtitle('Null probability vs covariate') +
				theme_bw() +
				theme(legend.position = "none")
		
		plot.list[['pval']] <- ggplot(data, aes(x = pi0.var, y = pvals, col = rejection)) + 
				geom_point(alpha = 0.75, size = 0.2) +
				xlab(covariate.name) +
				ylab(ylab.pval) +
				ggtitle(paste0('Target FWER level (alpha = ', alpha, ')')) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				theme_bw()
	} else if (covariate.type== 'categorical') {
		plot.list[['pi0']] <- ggplot(data, aes(x = pi0.var, y = pi0, col = rejection)) + 
				geom_boxplot(alpha = 0.75) +
				xlab(covariate.name) +
				ylab(ylab.pi0) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				ggtitle('Null probability vs covariate') +
				theme_bw() +
				theme(legend.position = "none")
		
		plot.list[['pval']] <- ggplot(data, aes(x = pi0.var, y = pvals, col = rejection)) + 
				geom_jitter(alpha = 0.75, size = 0.2, position = position_jitter(height = 0)) +
				xlab(covariate.name) +
				ylab(ylab.pval) +
				ggtitle(paste0('Target FWER level (alpha = ', alpha, ')')) +
				scale_colour_manual(values = c('darkgray', 'red')) +
				theme_bw()
	}
	
	obj <- plot_grid(plot.list[['pi0']], plot.list[['pval']], ncol = 2)
	
	cat('Generate combined PDF plot ...\n')
	pdf(paste0(file.name), width =10, height = 4)
	obj <- plot_grid(plot.list[['pi0']], plot.list[['pval']],  ncol=2)
	print(obj)
	dev.off()
	
	cat('Finished!\n')
	return(invisible(plot.list))
	
}
############################################################################################################
# Testing the informativeness of the covariate
# Copied from "DescTools" package to reduce dependency 
cochran.armitage.test <- function (x, alternative = c("two.sided", "increasing", "decreasing")) {
	DNAME <- deparse(substitute(x))
	if (!(any(dim(x) == 2))) 
		stop("Cochran-Armitage test for trend must be used with rx2-table", 
				call. = FALSE)
	if (dim(x)[2] != 2) 
		x <- t(x)
	nidot <- apply(x, 1, sum)
	n <- sum(nidot)
	Ri <- 1:dim(x)[1]
	Rbar <- sum(nidot * Ri)/n
	s2 <- sum(nidot * (Ri - Rbar)^2)
	pdot1 <- sum(x[, 1])/n
	z <- sum(x[, 1] * (Ri - Rbar)) / sqrt(pdot1 * (1 - pdot1) * 
					s2)
	STATISTIC <- z
	alternative <- match.arg(alternative)
	PVAL <- switch(alternative, two.sided = 2 * pnorm(abs(z), 
					lower.tail = FALSE), increasing = pnorm(z), decreasing = pnorm(z, 
					lower.tail = FALSE))
	PARAMETER <- dim(x)[1]
	names(STATISTIC) <- "Z"
	names(PARAMETER) <- "dim"
	METHOD <- "Cochran-Armitage test for trend"
	structure(list(statistic = STATISTIC, parameter = PARAMETER, 
					alternative = alternative, p.value = PVAL, method = METHOD, 
					data.name = DNAME), class = "htest")
}


#' A permutation test for assessing the informativeness of a numeric covariate for CAMT.
#'
#' The function implements a powerful statistical test for assessing the association between the p-value and a numeric covariate, exploiting the signal sparsity 
#' and potential nonlinearity. This is achieved by testing the association between two categorical variables after dichotomizing the p-values at the lower end 
#' and splitting the covariate into disjoint sets. An omnibus-type test is designed to combine evidence through various categorizations. Permutation is used
#' to assess the statistical significance. We recommend using the covariate if the p-value is highly significant (p < 0.005).
#'
#' @param pvals a numeric vector of p-values.
#' @param covariate a numeric vector of the covariate.
#' @param cutoffs a numeric vector of the cutoff points for dichotomizing the p-values. 
#' @param grps a vector of integer numbers specifying the numbers of breakpoints to divide the range of the covariate.
#' @param perm.no the number of permutation to assess the significance. Deafult is 999.
#' @param n.max an integer number specifying the maximum number of data points to be included. If the number of data points is larger than \code{n.max}, 
#' subsampling will be performed.
#' @param silence a logical value indicating whether to print out the process of the computation. 
#' 
#' @return A list with the elements
#' \item{stat.o}{the observed test statistic.}
#' \item{stat.p}{a vector of the test statistic from permutation.}
#' \item{p.value}{the omnibu p-value for testing the association.}
#' \item{x.cut.optim}{the optimal number of breakpoints for cutting the range of the covariate.}
#' \item{p.cut.optim}{the optimal cutoff to dichotomize the p-value.}

#' 
#' @author Jun Chen
#' @references Huang JY, ..., Zhang X, Chen J (2020) Leveraging biological and statistical covariates improves the detection power in epigenome-wide
#' association testing.21(88).
#' @keywords multiple testing, association test
#' @importFrom matrixStats rowMaxs colMaxs
#' @importFrom stats dgamma rgamma runif chisq.test quantile qchisq pchisq median rt
#' @examples
#'
#' data <- simulate.data(feature.no = 1000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Low', sig.strength = 'L3', cor.struct = 'None')
#' obj <- suppressWarnings(cov.test.n(data$pvals, data$pi0.var, perm.no = 2))
#' obj$p.value
#' 

#' 
#' @rdname cov.test.n
#' @export

cov.test.n <- function (pvals, covariate, cutoffs = quantile(pvals, c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)),
		grps = c(2, 4, 8, 16, 32), perm.no = 999, n.max = 100000, silence = TRUE) {
	
	save.seed <- get(".Random.seed", .GlobalEnv)
	
	n <- length(pvals)
	if (n > n.max) {
		ind <- sample(1:n, n.max)
		pvals <- pvals[ind]
		covariate <- covariate[ind]
	}
	
	
	stat.p.a <- stat.p.a1 <- stat.p.a2 <- array(NA, c(length(cutoffs), length(grps), perm.no))
	stat.o.a <- stat.o.a1 <- stat.o.a2 <- array(NA, c(length(cutoffs), length(grps)))
	
	for (i in 1:length(cutoffs)) {
		cutoff <- cutoffs[i]
		x <- as.numeric(pvals <= cutoff) 
		for (j in 1:length(grps)) {
			if (!silence) cat('.')
			grp <- grps[j]
			y <- cut(covariate, c(min(covariate) - 0.01, quantile(covariate, 1 / grp * (1:(grp-1))), max(covariate) + 0.01))
			
			if (!silence) cat('.')
			mat <- table(x, y)
			test.obj1 <- chisq.test(mat)
			test.obj2 <- cochran.armitage.test(mat)
			stat.o.a1[i, j] <- -pchisq(test.obj1$statistic, df = test.obj1$parameter, lower.tail = FALSE, log.p = TRUE)
			stat.o.a2[i, j] <- -log(2) - pnorm(abs(test.obj2$statistic), lower.tail = FALSE, log.p = TRUE)
			stat.o.a[i, j] <- max(c(stat.o.a1[i, j], stat.o.a2[i, j]))
			
			assign(".Random.seed", save.seed, .GlobalEnv)
			
			lpvs <- sapply(1:perm.no, function (k) {
						xp <- sample(x)		
						mat <- table(xp, y)
						test.obj1 <- chisq.test(mat)
						test.obj2 <- cochran.armitage.test(mat)
						c(-pchisq(test.obj1$statistic, df = test.obj1$parameter, lower.tail = FALSE, log.p = TRUE),
								-log(2) - pnorm(abs(test.obj2$statistic), lower.tail = FALSE, log.p = TRUE))
					})
			
			stat.p.a1[i, j, ]  <- lpvs[1, ]
			stat.p.a2[i, j, ]  <- lpvs[2, ]
			stat.p.a[i, j, ] <- pmax(lpvs[1, ], lpvs[2, ])
		}
	}
	
	stat.o1 <- colMaxs(stat.o.a)
	stat.o2 <- rowMaxs(stat.o.a)
	stat.o <- max(stat.o1)
	x.cut.optim <- grps[which.max(stat.o1)]
	p.cut.optim <- cutoffs[which.max(stat.o2)]
	stat.p <- apply(stat.p.a, 3, max)
	p.value <- mean(c(stat.p, stat.o) >= stat.o)
	
	stat.o1 <- max(stat.o.a1)
	stat.p1 <- apply(stat.p.a1, 3, max)
	p.value1 <- mean(c(stat.p1, stat.o1) >= stat.o1)
	
	stat.o2 <- max(stat.o.a2)
	stat.p2 <- apply(stat.p.a2, 3, max)
	p.value2 <- mean(c(stat.p2, stat.o2) >= stat.o2)
	
	return(list(stat.o = stat.o,  stat.p = stat.p, p.value = p.value,
#					stat.o1 = stat.o1,  stat.p1 = stat.p1, p.value1 = p.value1, 
#					stat.o2 = stat.o2,  stat.p2 = stat.p2, p.value2 = p.value2,
					x.cut.optim = x.cut.optim,  p.cut.optim = p.cut.optim))
}

#' A permutation test for assessing the informativeness of a categorical covariate for CAMT.
#'
#' The function implements a powerful statistical test for assessing the association between the p-value and a categorical covariate, exploiting the signal sparsity. 
#' This is achieved by testing the association between two categorical variables after dichotomizing the p-values at the lower end. An omnibus-type
#' test is designed to combine evidence through various dichotomization. Permutation is used to assess the statistical significance.   We recommend using the 
#' covariate if the p-value is highly significant (p < 0.005).
#'
#' @param pvals a numeric vector of p-values.
#' @param covariate a factor of the categorical covariate.
#' @param cutoffs a numeric vector of the cutoff points for dichotomizing the p-values. 
#' @param perm.no the number of permutation to assess the significance. Deafult is 999.
#' @param n.max an integer number specifying the maximum number of data points to be included. If the number of data points is larger than \code{n.max}, 
#' subsampling will be performed.
#' @param silence a logical value indicating whether to print out the process of the computation. 
#' 
#' @return A list with the elements
#' \item{stat.o}{the observed test statistic.}
#' \item{stat.p}{a vector of the test statistic from permutation.}
#' \item{p.value}{the omnibu p-value for testing the association.}
#' \item{p.cut.optim}{the optimal cutoff to dichotomize the p-value.}

#' 
#' @references Huang JY, ..., Zhang X, Chen J (2020) Leveraging biological and statistical covariates improves the detection power in epigenome-wide
#' association testing.21(88).
#' @keywords multiple testing, association test
#' @importFrom stats dgamma rgamma runif chisq.test quantile qchisq pchisq median rt
#' @examples
#'
#' data <- simulate.data(feature.no = 1000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Low', sig.strength = 'L3', cor.struct = 'None')
#' X <- factor(data$pi0.var >= median(data$pi0.var))
#' obj <- suppressWarnings(cov.test.c(data$pvals, X, perm.no = 2)) 
#' obj$p.value
#' 
#' @rdname cov.test.c
#' @export

cov.test.c <- function (pvals, covariate, cutoffs = quantile(pvals, c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)), 
		perm.no = 999, n.max = 100000, silence = TRUE) {
	
	save.seed <- get(".Random.seed", .GlobalEnv)
	
	n <- length(pvals)
	if (n > n.max) {
		ind <- sample(1:n, n.max)
		pvals <- pvals[ind]
		covariate <- covariate[ind]
	}
	
	y <- factor(covariate)
	
	stat.p.a <- array(NA, c(length(cutoffs), perm.no))
	stat.o.a <- array(NA, c(length(cutoffs)))
	
	for (i in 1:length(cutoffs)) {
		
		assign(".Random.seed", save.seed, .GlobalEnv)
		
		cutoff <- cutoffs[i]
		x <- as.numeric(pvals <= cutoff) 
		
		if (!silence) cat('.')
		test.obj <- chisq.test(table(x, y))
		stat.o.a[i] <- -pchisq(test.obj$statistic, df = test.obj$parameter, lower.tail = FALSE, log.p = TRUE)
		stat.p.a[i, ] <- sapply(1:perm.no, function (k) {
					xp <- sample(x)		
					test.obj.p <- chisq.test(table(xp, y))
					-pchisq(test.obj.p$statistic, df = test.obj.p$parameter, lower.tail = FALSE, log.p = TRUE)
				})
		
	}
	
	
	stat.o <- max(stat.o.a)
	p.cut.optim <- cutoffs[which.max(stat.o.a)]
	stat.p <- apply(stat.p.a, 2, max)
	p.value <- mean(c(stat.p, stat.o) >= stat.o)
	return(list(stat.o = stat.o, stat.p = stat.p, p.value = p.value, p.cut.optim = p.cut.optim))
}
############################################################################################################
# Simulation
Tmat <- function (feature.no = 10000, nb = 100, rho = 0) {
	# Simulate special correlation structure
	bs <- feature.no / nb
	mat <- diag(bs)
	mat[, ] <- rho
	diag(mat) <- 1
	obj <- eigen(mat)
	T1 <- obj$vectors %*% diag(sqrt(obj$values)) %*% t(obj$vectors)
	
	mat <- diag(bs)
	mat[, ] <- -rho
	mat[1:(bs/2), 1:(bs/2)] <- rho
	mat[(bs/2+1):bs, (bs/2+1):bs] <- rho
	diag(mat) <- 1
	obj <- eigen(mat)
	T2 <- obj$vectors %*% diag(sqrt(obj$values)) %*% t(obj$vectors)
	
	return(list(T1=T1, T2=T2))
}

paras.mapping.func <- function (
		sig.densities = c(Low=3.5, Medium=2.5, High=1.5),
		sig.strengths = log(seq(exp(2), exp(3), len=8)),
		pi.strengths = c(None=0.0, Moderate=1.0, Strong=1.5),
		f1.strengths = c(None=0.0, Moderate=0.25, Strong=0.5)
) {
	
	names(sig.strengths) <- paste0('L', 1:8)
	names(sig.densities) <- c('Low', 'Medium', 'High')
	names(pi.strengths) <- c('None', 'Moderate', 'Strong')
	names(f1.strengths) <- c('None', 'Moderate', 'Strong')
	return(list(sig.strengths = sig.strengths, sig.densities = sig.densities,
					pi.strengths = pi.strengths, f1.strengths = f1.strengths))
	
}

#' Simulate p-values and covariates under various scenarios.
#'
#' The function simulates p-values and covariates under different signal structures (density and strength) and covariate models. 
#' It also allows special correlation structures among the p-values.
#'
#' @param paras.mapping a list with slots \code{'sig.strengths', 'sig.densities', 'pi.strengths', 'f1.strengths'}, which define the actual
#' numeric values used for different levels of \code{'sig.strength', 'sig.density', 'covariate.strength'}.  \code{'covariate.strength'} consists
#' of both \code{'pi.strengths'} and \code{'f1.strengths'}. The mapping can be generated by using the \code{paras.mapping.func} function.
#' @param covariate.strength a character string from \code{'None', 'Moderate', 'Strong'} indicating the covariate strength.
#' @param covariate.model a character string from \code{'pi0', 'f1', 'both'} indicating whether the prior null proability, 
#' the alterantive distribution or both are affected by the covariate.
#' @param covariate.dist a character string from \code{'Normal', 'Uniform', 'T'} indicating the distribution of the covariate.
#' For t-distribution, a degree of freedom of five is used.
#' @param null.model a character string from \code{'Unif', 'Left', 'Right'} indicating whether the null distribution of the p-value is uniform, 
#' left skewed or right skewed.
#' @param skewness a numeric value indicating the skewness of the p-value distribution under the null (mean of the z-value)
#' @param f1.sd a numeric value indicating the variance of the z-score under the alternative. Default is 1. 
#' @param feature.no an integar, the number of features to be simulated.
#' @param sig.dist  a character string from \code{'Normal', 'Gamma'} indicating the distribution of the z-value under alternative.
#' @param sig.density  a character string from \code{'Low', 'Medium', 'High'} indicating the level of signal density.
#' @param sig.strength a character string from \code{'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'} indicating the level of signal strength.
#' @param cor.struct a character string from \code{'None', 'Block', 'AR1'} indicating the correlation structure to be simulated including
#' no correlation, block correlation and AR(1) correlation. 
#' @param cor.rho  a numeric value giving the correlation coefficient for \code{'Block', 'AR1'}.
#' @param cor.sign  a character string from \code{'PosCor', 'PosNegCor'} indicating whether to simulate only positive or both postive and negative correlations. 
#' @param cor.nblock an integar, the number of blocks to be simulated for block correlation structure
#'
#' @return A list with the elements
#' \item{pvals}{a numeric vector of p-values.}
#' \item{pi0.var}{a vector of covariate values for the prior null probability.}
#' \item{f1.var}{a vector of covariate values for the alternative distribution.}
#' \item{pi0}{ a vector of the simulated null probabilities.}
#' \item{truth}{a vector simulated truth for H0 (=0) or H1 (=1).}
#' \item{lfdr}{a vector of oracle LFDR based on the simulated pi0, f1.}
#' \item{fdr}{a vector of oracle FDR based on LFDR.}
#' 
#' @author Jun Chen
#' @references Xianyang Zhang, Jun Chen. Covariate Adaptive False Discovery Rate Control with Applications to Omics-Wide Multiple Testing. JASA. To appear.
#' @keywords FDR, multiple testing
#' @importFrom stats rnorm dnorm pnorm arima.sim rbinom
#' @importFrom cowplot plot_grid
#' @examples
#'
#' data <- simulate.data(feature.no = 10000, covariate.strength = 'Moderate', covariate.model = 'pi0',
#'		sig.density = 'Medium', sig.strength = 'L3', cor.struct = 'None')
#' camt.fdr.obj <- camt.fdr(pvals = data$pvals, pi0.var = data$pi0.var, f1.var = data$f1.var, 
#'		alg.type = 'EM', control.method = 'knockoff+')
#' plot.camt.fdr(camt.fdr.obj, covariate = as.vector(rank(data$pi0.var)), covariate.name = 'Covariate rank',
#'		log = TRUE, file.name = 'CovariateModerate.pdf')
#' 
#' @rdname simulate.data
#' @export
simulate.data <- function (
		paras.mapping       = paras.mapping.func(),
		covariate.strength  = c('None', 'Moderate', 'Strong'),  
		covariate.model     = c('pi0', 'f1', 'both'),
		covariate.dist      = c('Normal', 'Uniform', 'T'),
		null.model          = c('Unif', 'Left', 'Right'),
		skewness            = 0.15,
		f1.sd               = 1.0,
		feature.no          = 10000,
		sig.dist            = c('Normal', 'Gamma'), 
		sig.density         = c('Low', 'Medium', 'High'),
		sig.strength        = c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'),
		cor.struct          = c('None', 'Block', 'AR1'),
		cor.rho             = 0,
		cor.sign            = c('PosCor', 'PosNegCor'),
		cor.nblock          = 500) {

	
	sig.dist <- match.arg(sig.dist)
	sig.density <- match.arg(sig.density)
	sig.strength <- match.arg(sig.strength)
	covariate.strength <- match.arg(covariate.strength)
	covariate.model <- match.arg(covariate.model)
	covariate.dist <- match.arg(covariate.dist)
	null.model <- match.arg(null.model)
	cor.struct <- match.arg(cor.struct)
	cor.sign <- match.arg(cor.sign)
	
	# These could be altered for other simulation purpose

	sig.densities <- paras.mapping[['sig.densities']]
	sig.strengths <- paras.mapping[['sig.strengths']]
	pi.strengths <- paras.mapping[['pi.strengths']]
	f1.strengths <- paras.mapping[['f1.strengths']]
	
	if (covariate.model %in% c('pi0')) {
		pi.k <- pi.strengths[covariate.strength]
		f1.k <- 0
	} 
	if (covariate.model %in% c('f1')) {
		pi.k <- 0
		f1.k <- f1.strengths[covariate.strength]
	} 
	if (covariate.model %in% c('both')) {
		pi.k <- pi.strengths[covariate.strength]
		f1.k <- f1.strengths[covariate.strength]
	} 
	
	# Generate covariate for pi0 and f1
    # We simulate different covariates
	if (covariate.dist == 'Normal') {
		x1 <- rnorm(feature.no)
		x2 <- rnorm(feature.no)
	} 
	
	if (covariate.dist == 'Uniform') {
		x1 <- scale(runif(feature.no))
		x2 <- scale(runif(feature.no))
	}
	
	if (covariate.dist == 'T') {
		x1 <- rt(feature.no, df = 5)
		x2 <- rt(feature.no, df = 5)
	}
	
	eta <- sig.densities[sig.density] + x1 * pi.k
	pi <- exp(eta) / (1 + exp(eta))
	
	truth <- rbinom(feature.no, prob = 1 - pi, size = 1)
	ns <- sum(truth)
	
	# Generate the correlated signals under the null
	if (cor.struct != 'None') {
		if (null.model != 'Unif') {
			cat('Correlated signals will be simulated under the uniform null only! "null.model" will be automatically reset to be "Unif". \n')
		} else {
			if (cor.struct == 'Block') {
				obj <- Tmat(feature.no, cor.nblock, cor.rho)
				bs <- feature.no / cor.nblock
				if (cor.sign == 'PosCor') {
					stats <- as.vector(obj$T1 %*% matrix(rnorm(feature.no), bs, cor.nblock))
					
				}
				if (cor.sign == 'PosNegCor') {
					stats <- as.vector(obj$T2 %*% matrix(rnorm(feature.no), bs, cor.nblock))
				}
			} 
			if (cor.struct == 'AR1') {
				if (cor.sign == 'PosCor') {
					stats <- as.vector(arima.sim(n = feature.no, list(ar = c(cor.rho)), n.start = 100000) * sqrt(1 - cor.rho^2))
				}
				if (cor.sign == 'PosNegCor') {
					stats <- as.vector(arima.sim(n = feature.no, list(ar = c(-cor.rho)), n.start = 100000) * sqrt(1 - cor.rho^2))
				}
			}
			skewness <- 0
		}
		
	} else {
		if (null.model == 'Unif') {
			stats <- rnorm(feature.no)
			skewness <- 0
		}
		if (null.model == 'Left') {
			stats <- rnorm(feature.no, -abs(skewness), 1)
			skewness <-  -abs(skewness)
		}
		if (null.model == 'Right') {
			stats <- rnorm(feature.no, abs(skewness), 1)
			skewness <- abs(skewness)
		}
		
	}
	
	if (sig.dist == 'Normal') {
		f1.mean <- sig.strengths[sig.strength] * (exp(f1.k * x2) / (1 + exp(f1.k * x2))) * 2
				
		if (null.model == 'Unif') {
			stats[truth == 1] <- stats[truth == 1] * f1.sd + f1.mean[truth == 1]
		} else {
			stats[truth == 1] <- rnorm(ns) * f1.sd + f1.mean[truth == 1]
		}
		
		lfdr <- pi * dnorm(stats, mean = skewness) / (pi * dnorm(stats, mean = skewness) + (1 - pi) * dnorm(stats, f1.mean, sd = f1.sd))
		
		out <- sort(lfdr, index.return = TRUE)
		fdr <- cumsum(out$x) / (1:(length(lfdr)))
		fdr <- fdr[order(out$ix)]
	}
	
	if (sig.dist == 'Gamma') {
		
		# Under gamma distribution (a) the nulls and non-nulls are not correlated (b) we do not shrink the variance in the covariate model 'f1' and 'both'
		f1.mean <- sig.strengths[sig.strength] * (exp(f1.k * x2) / (1 + exp(f1.k * x2))) * 2
		
		# Match the mean and variance of the normal counterpart
		stats[truth == 1] <- rgamma(ns, shape = 2, scale = 1 / sqrt(2)) - sqrt(2) + f1.mean[truth == 1]
		
		# Calculate local false discovery rate
		lfdr <- pi * dnorm(stats, mean = skewness) / (pi * dnorm(stats, mean = skewness) + (1 - pi) * 
					dgamma(stats - f1.mean + sqrt(2), shape = 2, scale = 1 / sqrt(2))) 
		
		out <- sort(lfdr, index.return = TRUE)
		fdr <- cumsum(out$x) / (1:(length(lfdr)))
		fdr <- fdr[order(out$ix)]
		
	}
	
	# Generate p values - one sided p-value
	pvals <- (1 - pnorm(stats))
	
	return(list(pvals = pvals, pi0.var = x1, f1.var = x2, pi0 = pi, truth = truth,  fdr = fdr, lfdr = lfdr))
}

