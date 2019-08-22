## 
Solve_q_ordered_simple <- function(pvals, tau, eps, ADMM_params){
	target.num <- 5000
	n <- length(pvals)
	if (n <= target.num){
		qhat <- Solve_q_ordered(pvals, tau, eps, ADMM_params) 
		return(qhat)
	}
	num.reps <- ceiling(n / target.num)
	new.pvals <- sapply(1:target.num, function(i){
				ind <- pmin((i - 1) * num.reps + (1:num.reps), n)
				pvals[ind[1]]
			})
	qhat <- Solve_q_ordered(new.pvals, tau, eps, ADMM_params)
	qhat <- rep(qhat, each = num.reps)[1:n]
	return(qhat)
}


runtime.calc <- function (methods = 'CAMT', ns = 1000, iters = 1, fdr.cutoff = 0.05) {
	res <- array(NA, c(length(methods), length(ns), length(iters)), dimnames=list(methods, paste(ns), paste(iters)))
	for (n in ns) {
		cat('*')
		for (method in methods) {
			cat('!')
			for (i in iters) {
				cat('.')
				data <- simulate.data(feature.no = n, covariate.strength = 'Moderate', sig.density = 'Medium', sig.strength = 'L3')
				
				pvalue <- data$pvals
				x1 <- data$pi0.var
				x2 <- data$f1.var
				
				if (method == 'CAMT') {
					runtime <- system.time({
								camt.fdr(pvals = pvalue, pi0.var = x1, f1.var = x2, 
										alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)
							})[3]
				}
				
				if (method == 'IHW') {
					runtime <- system.time({
								ihw.res <- ihw(pvalue ~ x1, alpha = fdr.cutoff)
							})[3]
				}
				
				
				if (method == 'FDRreg') {
					runtime <- system.time({
								FDRreg.res <- FDRreg(CAMT:::ztransform(pvalue), matrix(x1), nulltype = 'theoretical')
							})[3]
				}
				
				
				if (method == 'SABHA') {
					runtime <- system.time({
								index0 <- order(x1)
								index1 <- order(index0)
								pvals <- pvalue[index0]
								qhat <- Solve_q_step(pvals, 0.5, 0.1)
								SABHA.res <- SABHA_method(pvals, qhat, fdr.cutoff, 0.5)	
							})[3]
					
				}
				
				ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr
				
				
				if (method == 'SABHA2') {
					runtime <- system.time({
								index0 <- order(x1)
								index1 <- order(index0)
								pvals <- pvalue[index0]
								qhat_ordered <- Solve_q_ordered_simple(pvals, 0.5, 0.1, ADMM_params)
								SABHA.res <- SABHA_method(pvals, qhat_ordered, fdr.cutoff, 0.5)	
							})[3]
					
				}
				
				if (method == 'AdaPT') {
					runtime <- system.time({
								x.df <- data.frame(x1 = x1, x2 = x2)
								pi_formulas <- paste0("x1")
								mu_formulas <- paste0("x2")
								adapt1.res <- adapt_glm(x = x.df, pvals = pvalue, pi_formulas = pi_formulas, mu_formulas = mu_formulas,
										verbose = list(print = FALSE, fit = FALSE, ms = FALSE))
							})[3]
					
				}
				
				if (method == 'BL') {
					runtime <- system.time({
								pi0x <- lm_pi0(pValues = pvalue, X = x1, smooth.df = 1)
								BL.res <- p.adjust(pvalue, method = "fdr") * pi0x$pi0
							})[3]
				}
				res[method, paste(n), paste(i)] <- runtime
				write(runtime, file=paste0(method, '.', n, '.', i, '.txt'))
				
			}
		}
	}
}


