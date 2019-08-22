
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
			dimnames=list(struct.types, prior.strengths, paste(feature.nos), sig.dists, sig.densities, paste(sig.strengths), methods, measures))	
	
	for (struct.type in struct.types) {
		for (prior.strength in prior.strengths) {
			cat('$')
			for (feature.no in feature.nos) {
				cat('*')
				for (sig.dist in sig.dists) {
					cat('!')
					for (sig.density in sig.densities) {
						cat('%')

						# Generate data
						data <- simulate.data(paras.mapping = CAMT:::paras.mapping.func(sig.densities = c(Low=3.5, Medium=2.5, High=1.5)),
								covariate.strength = prior.strength, sig.density=sig.density, sig.strength = 'L1',
								sig.dist = sig.dist, feature.no = as.numeric(feature.no), 
								null.model = null.model, covariate.model = covariate.model, covariate.dist = covariate.dist,
								cor.struct = cor.struct, cor.nblock = cor.nblock, cor.sign = cor.sign, cor.rho = cor.rho)
						
						pvalue <- data$pvals
						pvalue[pvalue <= 1e-15] <- 1e-15
						truth <- data$truth
						x1 <- data$pi0.var
						x2 <- data$f1.var
						
						res.obj <- list()
						
						for (j in 1:length(sig.strengths))  {
							cat('#')
							# Generate data
							fdr.cutoff <- sig.strength <- sig.strengths[j]
							
							
							if (struct.type == 'Covariate') {
								cat('.')
								
								res.obj[['Oracle']] <- list(fdr=data$fdr)
								
								try({
											if (j == 1) {
												res.obj[['CAMT']] <- list(fdr = camt.fdr(pvals = pvalue, pi0.var = x1, f1.var = x1, 
																alg.type = 'EM', control.method = 'knockoff+', trace = FALSE, pi0.low = 0.1)$fdr)
											}
											cat('.')
											
											
											
										})
								try({
											
											if (j == 1) {
												cat('.')
												ihwRes <- ihw(pvalue ~ x1, alpha = fdr.cutoff)
												res.obj[['IHW']] <- list(fdr=adj_pvalues(ihwRes))
											}
											
											
										})
								
								try({
											if (j == 1) {
												cat('.')
												res.obj[['FDRregT']] <- list(fdr=FDRreg(CAMT:::ztransform(pvalue), matrix(x1), nulltype = 'theoretical')$FDR)
											}
											
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
											if (j == 1) {
												cat('.')
												index0 <- order(x1)
												index1 <- order(index0)
												
												fdr <- rep(1, feature.no)
												pvals <- pvalue[index0]
												qhat <- Solve_q_step(pvals, 0.5, 0.1)
											}
											
											
											#	ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3)
											#   qhat <- Solve_q_ordered(pvals, tau=0.5, eps=0.1, ADMM_params)
											
											
											
											fdr[SABHA_method(pvals, qhat, fdr.cutoff, 0.5)]  <- fdr.cutoff / 2
											res.obj[['SABHA']] <- list(fdr=fdr[index1])
										})
								
							}
							
							try({
										if (j == 1) {
											cat('.')
											pi0x <- lm_pi0(pValues = pvalue, X = x1, smooth.df = 1)
											res.obj[['BL']] <- list(fdr = p.adjust(pvalue, method = "fdr") * pi0x$pi0)	
										}
										
									})
							
							
							if (j == 1) {
								cat('.')
								res.obj[['BH']] <- list(fdr=p.adjust(pvalue, 'fdr'))
							}
							
							if (j == 1) {
								cat('.')
								res.obj[['ST']] <- list(fdr=qvalue(pvalue)$qvalue)
							}
							
							
							for (method in methods) {
								
								# Evaluation
								if (is.null(res.obj[[method]])) {
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, paste(sig.strength), method, c('FDR', 'Power')] <- NA
									next
								} 
								
								fdr <- res.obj[[method]]$fdr
								
								if (sum(fdr <= fdr.cutoff) == 0) {
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, paste(sig.strength), method, 'FDR'] <- 0
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, paste(sig.strength), method, 'Power'] <- 0
								} else {
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, paste(sig.strength), method, 'FDR'] <- mean(truth[fdr <= fdr.cutoff] == 0)
									res[struct.type, prior.strength, paste(feature.no), sig.dist, sig.density, paste(sig.strength), method, 'Power'] <- sum(truth[fdr <= fdr.cutoff])
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

prefix <- 'cov_norm_k_0_r_0.v2.100000'
resdir <-  file.path("~/project/covfdr/", prefix)
source(file.path("~/project/covfdr/", prefix, "Cluster_mayo.R"))

# Covariate Structure
paras <- list()
paras$resdir <- resdir
paras$prefix <- prefix
paras$fdr.cutoff <- 0.05
paras$feature.nos <- '100000'
paras$struct.types <- c('Covariate')
paras$sig.dists <- c('Normal')
paras$prior.strengths <- c('Moderate')
paras$sig.densities <- c('Low')
paras$sig.strengths <- c(0.01, 0.05, 0.10, 0.15, 0.20)
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
