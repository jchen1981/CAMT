
func <- function(part, paras) {
	require(qvalue)
	require(IHW)
	require(FDRreg)
	require(adaptMT)
	require(CAMT)
	require(splines)
	require(swfdr)
	source('https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R')
	
	set.seed(part)
	feature.nos <- paras$feature.nos
	sig.dists <- paras$sig.dists
	struct.types <- paras$struct.types
	prior.strengths <- paras$prior.strengths
	sig.densities <- paras$sig.densities
	sig.strengths <- paras$sig.strengths
	fdr.cutoff <- paras$fdr.cutoff
	null.model <- paras$null.model
	covariate.model <- paras$covariate.model
	covariate.dist <- paras$covariate.dist
	
	cor.struct <- paras$cor.struct
	cor.nblock <- paras$cor.nblock 
	cor.rho <- paras$cor.rho 
	cor.sign <- paras$cor.sign 
	
	resdir <- paras$resdir
	prefix <- paras$prefix
	
	resdir <- gsub('/$', '', resdir)
	
	sink(file.path(resdir, paste(prefix, "_",  part, ".log", sep="")))
	cat(date(), '\n')
	
	methods <- c('Oracle',  'CAMT',  'IHW', 'FDRregT', 'FDRregE', 'AdaPT', 'SABHA', 'BL',  'BH', 'ST')

	measures <- c('FDR', 'Power')
	res <- array(NA, c(length(struct.types), length(prior.strengths), length(feature.nos), length(sig.dists), length(sig.densities), length(sig.strengths),
					length(methods), length(measures)), 
			dimnames=list(struct.types, prior.strengths, feature.nos, sig.dists, sig.densities, sig.strengths, methods, measures))	
	
	for (struct.type in struct.types) {
		for (prior.strength in prior.strengths) {
			cat('$')
			for (feature.no in feature.nos) {
				cat('*')
				for (sig.dist in sig.dists) {
					cat('!')
					for (sig.density in sig.densities) {
						cat('%')
						for (sig.strength in sig.strengths) {
							cat('#')
							# Generate data
							data <- simulate.data(covariate.strength = prior.strength, sig.density=sig.density, sig.strength = sig.strength,
									sig.dist = sig.dist, feature.no = feature.no, 
									null.model = null.model, covariate.model = covariate.model, covariate.dist = covariate.dist,
									cor.struct = cor.struct, cor.nblock = cor.nblock, cor.sign = cor.sign, cor.rho = cor.rho)
							
							pvalue <- data$pvals
							pvalue[pvalue <= 1e-15] <- 1e-15
							truth <- data$truth
							feature.no <- length(pvalue)
							x1 <- data$pi0.var
							x2 <- data$f1.var
							
							res.obj <- list()

							
							if (struct.type == 'Covariate') {
								cat('.')
								
								res.obj[['Oracle']] <- list(fdr=data$fdr)
								
								try({
											cat('.')
											res.obj[['CAMT']] <- list(fdr = camt.fdr(pvals = pvalue, pi0.var = x1, f1.var = x1, 
															alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)$fdr)
											
										})
								try({
											cat('.')
											ihwRes <- ihw(pvalue ~ x1, alpha = fdr.cutoff)
											res.obj[['IHW']] <- list(fdr=adj_pvalues(ihwRes))
											
										})
								
								try({
											cat('.')
											res.obj[['FDRregT']] <- list(fdr=FDRreg(CAMT:::ztransform(pvalue), matrix(x1), nulltype = 'theoretical')$FDR)
										})	
								
								try({
											cat('.')
											res.obj[['FDRregE']] <- list(fdr=FDRreg(CAMT:::ztransform(pvalue), matrix(x1), nulltype = 'empirical')$FDR)
										})	
								
								try({
											cat('.')
											fdr <- rep(1, feature.no)
											x.df <- data.frame(x = x1)
											formulas <- paste0("x")
											adapt1.res <- adapt_glm(x = x.df, pvals = pvalue, pi_formulas = formulas, mu_formulas = formulas, verbose = list(print = FALSE, fit =
																	FALSE, ms = FALSE))
											fdr[pvalue <= adapt1.res$s[, floor(fdr.cutoff * 100)]]  <- fdr.cutoff / 2
											res.obj[['AdaPT']] <- list(fdr=fdr)
											
										})
								try({ 
											cat('.')
											index0 <- order(x1)
											index1 <- order(index0)
											
											fdr <- rep(1, feature.no)
											pvals <- pvalue[index0]
											
											#	ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3)
			                                #   qhat <- Solve_q_ordered(pvals, tau=0.5, eps=0.1, ADMM_params)
			
											qhat <- Solve_q_step(pvals, 0.5, 0.1)
											
											fdr[SABHA_method(pvals, qhat, fdr.cutoff, 0.5)]  <- fdr.cutoff / 2
											res.obj[['SABHA']] <- list(fdr=fdr[index1])
										})
								
							}
							
							try({
										cat('.')
										pi0x <- lm_pi0(pValues = pvalue, X = x1, smooth.df = 1)
										res.obj[['BL']] <- list(fdr = p.adjust(pvalue, method = "fdr") * pi0x$pi0)	
									})

							
							cat('.')
							res.obj[['BH']] <- list(fdr=p.adjust(pvalue, 'fdr'))
							
							cat('.')
							res.obj[['ST']] <- list(fdr=qvalue(pvalue)$qvalue)
							
							for (method in methods) {
								
								# Evaluation
								if (is.null(res.obj[[method]])) {
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, c('FDR', 'Power')] <- NA
									next
								} 
								
								fdr <- res.obj[[method]]$fdr
								
								if (sum(fdr <= fdr.cutoff) == 0) {
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'FDR'] <- 0
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'Power'] <- 0
								} else {
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'FDR'] <- mean(truth[fdr <= fdr.cutoff] == 0)
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, sig.strength, method, 'Power'] <- sum(truth[fdr <= fdr.cutoff]) / sum(truth)
								}
							}
							
						}					
					}
				}
			}
		}
	}
	
	warnings()
	cat('\n', date(), '\n')
	sink()
	save(res, file=file.path(resdir, paste(prefix, "_res",  part, ".Rdata", sep="")))
	return(res)
}

prefix <- 'cov_norm_k_0_r_0_1.v2'
resdir <-  file.path("~/project/covfdr/", prefix)
source(file.path("~/project/covfdr/", prefix, "Cluster_mayo.R"))

# Covariate Structure
paras <- list()
paras$resdir <- resdir
paras$prefix <- prefix
paras$fdr.cutoff <- 0.05
paras$feature.nos <- 10000
paras$struct.types <- c('Covariate')
paras$sig.dists <- c('Normal')
paras$prior.strengths <- c('None', 'Moderate', 'Strong')
paras$sig.densities <- c('Low', 'Medium', 'High')
paras$sig.strengths <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')
paras$null.model <- 'Unif'
paras$covariate.model <- 'pi0'
paras$covariate.dist <- 'Normal'
paras$cor.struct <- 'None'
paras$cor.nblock <- 500
paras$cor.rho <- 0
paras$cor.sign <- 'PosCor'

setwd(resdir)
res <- clsapply(1:100, func, paras, queque='1-day', tempdir=file.path(resdir, 'tmpC'))
setwd(resdir)
save(res, file=file.path(resdir, paste(prefix, "_res.Rdata", sep="")))
