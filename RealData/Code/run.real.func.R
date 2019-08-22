run.all.method <- function (pvals, covariate, alpha.list = seq(0.01, 0.2, by = 0.01)) {
	
	x1 <- covariate
    x1.ns <- ns(x1, df=6)
	
	# Run algorithms
	cat('Runing CAMT ...\n')
	CAMT.obj <- camt.fdr(pvals = pvals, pi0.var = x1.ns, f1.var = x1.ns, 
			alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)
	
	cat('Runing IHW ...\n')
	IHW.obj <- adj_pvalues(ihw(pvals ~ x1, alpha = 0.05))
	
	cat('Runing FDRreg ...\n')
	FDRreg.obj <- FDRreg(CAMT:::ztransform(pvals), x1.ns, nulltype = 'theoretical')
	
	cat('Runing SABHA ...\n')
	qhat_step <- Solve_q_step(pvals, tau = 0.5, eps = 0.1)
	
	cat('Runing AdaPT...\n')
	x.df <- data.frame(x = x1)
	formulas <- paste0("ns(x, df=6)")
	AdaPT.obj <- adapt_glm(x = x.df, pvals = pvals, pi_formulas = formulas, mu_formulas = formulas, 
			verbose = list(print = FALSE, fit = FALSE, ms = FALSE))
	
	cat('Runing BL ...\n')
	pi0x <- lm_pi0(pValues = pvals, X = x1, smooth.df = 6)
	BL.obj <- p.adjust(pvals, method = "fdr") * pi0x$pi0
	
	cat('Runing BH/ST ...\n')
	ST.obj <- qvalue(pvals)
	BH.obj <- p.adjust(pvals, 'fdr')
	
	n <- length(pvals)
	num_alpha <- length(alpha.list)
	max_alpha <- max(alpha.list)
	
	methods <- c('CAMT',  'IHW', 'FDRreg', 'AdaPT', 'SABHA', 'BL', 'BH', 'ST')
	NumRej <- matrix(0, nrow = length(methods), ncol = num_alpha)
	rownames(NumRej) <- methods
	
	for(i in 1:num_alpha){
		NumRej['CAMT', i] <- sum(CAMT.obj$fdr <= alpha.list[i])
		NumRej['IHW', i] <- sum(IHW.obj <= alpha.list[i])
		NumRej['FDRreg', i] <- sum(FDRreg.obj$FDR <= alpha.list[i])
		NumRej['AdaPT', i] <- AdaPT.obj$nrejs[i]
		NumRej['SABHA', i] <- length(SABHA_method(pvals, qhat_step, alpha.list[i], tau=0.5))
		NumRej['BL', i] <- sum(BL.obj <= alpha.list[i])
		NumRej['BH', i] <- sum(BH.obj <= alpha.list[i])
		NumRej['ST', i] <- sum(ST.obj$qvalue <= alpha.list[i])
	}
	
	return(NumRej)

}


