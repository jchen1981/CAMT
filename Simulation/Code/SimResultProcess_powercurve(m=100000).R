# Set the working directory to that contains the simulation R data objects

require(reshape)
require(ggplot2)


for (prefix in c('cov_norm_k_0_r_0.v2.100000')) {	
	cat(prefix)
	load(paste0(prefix, '_res.Rdata'))
	struct.types <- c('Covariate')
	if (grepl('norm', prefix)){
		sig.dists <- c('Norm')
	} else {
		sig.dists <- c('Gamma')
	}
	###################################################################################################
# Combine the results
	feature.nos <- '100000'
# Sequences should be exactly the same as in the simulation
	prior.strengths <- c('Strong')
	sig.densities <- c('Low')
	sig.strengths <-  paste(c(0.01, 0.05, 0.10, 0.15, 0.20))
	methods <- c('Oracle',  'CAMT',  'IHW', 'FDRregT', 'FDRregE', 'AdaPT', 'SABHA', 'BL',  'BH', 'ST')
	measures <- c('FDR', 'Power')
	res.a <- array(NA, c(length(struct.types), length(prior.strengths), length(feature.nos), length(sig.dists), length(sig.densities), length(sig.strengths),
					length(methods), length(measures), length(res)), 
			dimnames=list(StructType=struct.types, OrderPrior=prior.strengths, FeatureNumber=feature.nos, SignalDist=sig.dists,
					SignalDensity=sig.densities, SignalStrength=paste(sig.strengths), Method=methods, 
					Measure=measures, Iteration=paste(1:length(res))))	
	
	for (i in 1:length(res)) {
		res.a[, , , , , , , , i] <- res[[i]]
	}
	
	###################################################################################################
# Calculate the standard error
	res.df <- melt(res.a)
	colnames(res.df)[ncol(res.df)] <- 'Value'
	
	m  <- aggregate(Value ~ StructType + OrderPrior + FeatureNumber + SignalDist + SignalDensity + SignalStrength + Method + 
					Measure, res.df, function(x) mean(x[!is.na(x)]))
	se <- aggregate(Value ~ StructType + OrderPrior + FeatureNumber + SignalDist + SignalDensity + SignalStrength + Method + 
					Measure, res.df, function(x) {
				ind <- !is.na(x)
				sd(x[ind]) / sqrt(length(x[ind]))})
	sd <- aggregate(Value ~ StructType + OrderPrior + FeatureNumber + SignalDist + SignalDensity + SignalStrength + Method + 
					Measure, res.df, function(x) {
				ind <- !is.na(x)
				sd(x[ind])})
	
	ymin <- m[, ncol(m)] - 1.96 * se[, ncol(m)]
	ymin <- ifelse(ymin > 0, ymin, 0)
	ymax <- m[, ncol(m)] + 1.96 * se[, ncol(m)]
	res.df2 <- cbind(m, SD = sd[, ncol(sd)], ymax=ymax, ymin=ymin,  SE=se[, ncol(se)], ErrRate=error.info[, ncol(m)])
	
	res.df2$FeatureNumber <- factor(res.df2$FeatureNumber, levels=feature.nos)
	res.df2$SignalDensity <- factor(res.df2$SignalDensity, levels=sig.densities)
	res.df2$SignalStrength <- factor(res.df2$SignalStrength, levels=sig.strengths)
	res.df2$OrderPrior <- factor(res.df2$OrderPrior, levels=prior.strengths)
	res.df2$Method <- factor(res.df2$Method, levels=c('Oracle', 'CAMT',  'AdaPT',  'IHW', 'FDRregT', 'FDRregE', 'SABHA', 'BL', 'ST',  'BH'))
	
	###################################################################################################
# Plot - Annonation, theme
	algs <- c('Oracle', 'CAMT',  'AdaPT',  'IHW', 'FDRregT', 'FDRregE', 'SABHA', 'BL', 'ST',  'BH')
	algs.ann <- c('Oracle', 'CAMT',  'AdaPT',  'IHW', 'FDRreg(T)', 'FDRreg(E)', 'SABHA', 'BL', 'ST',  'BH')
	names(algs.ann) <- algs
	
	nMeth <- length(algs)
	cols <- c('Oracle'="#999999",'CAMT'="red",  'AdaPT'= 'blue','IHW'="#D55E00", 'FDRregT'="#CC79A7", 
			'FDRregE'="#CC79A7", 'SABHA'="#009E73", 'BL' = "#0072B2", 'BH'="gray60", 'ST'="gray30")
	ltys <-  c('Oracle'=1, 'CAMT'=1, 'AdaPT'= 2,'IHW'= 1,  'FDRregT'=2, 'FDRregE'=2, 'SABHA'=1, 'BL'=2,'BH'=2, 'ST'=1)
	shapes <- c('Oracle'=17,'CAMT'=19,   'AdaPT'=3,  'IHW'=2, 'FDRregT'=0, 'FDRregE'=0,'SABHA'=4, 'BL'=6, 'BH'=15, 'ST'=5)
	
	names(cols) <- names(shapes) <- names(ltys) <- algs.ann[names(cols)]
	
	sig.density.ann <- c('Low'='Sparse Signal', 'Medium'='Medium Signal', 'High'='Dense Signal')
	prior.strength.ann <- c('None'='Non-informative Prior', 'Moderate' = 'Moderate Prior', 'Strong'='Strong Prior')
	
	###################################################################################################
# Methods to be plotted and their orders
	
	for (Struct in struct.types) {
		if (Struct == 'Covariate') {
			algs2 <-  c('Oracle',  'CAMT', 'AdaPT', 'IHW', 'FDRregT', 'SABHA',  'BL', 'BH', 'ST')
			algs2.ann <- algs.ann[algs2]
			cols2 <- cols[algs2.ann]
			ltys2 <- ltys[algs2.ann]
			shapes2 <- shapes[algs2.ann]
		}

		
		dodge <- position_dodge(width=0.9)
		for (Dist in sig.dists) {
			cat("*")
			plot.list <- list()
			
			for (Mea in c('FDR', 'Power')) {
				cat(".")
				
				res3 <- subset(res.df2, StructType %in% Struct & Measure %in% Mea &  SignalDist %in% Dist & Method %in% algs2, drop=TRUE)	
				
				levels(res3$Method) <- algs.ann[levels(res3$Method)]
				levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
				
				res3$SignalStrength <- as.numeric(as.character(res3$SignalStrength))
				levels(res3$OrderPrior) <- prior.strength.ann[levels(res3$OrderPrior)]

				obj <- ggplot(res3, aes(x=SignalStrength, y=Value,  group=Method, color=Method, linetype=Method, shape=Method)) 

				obj <- obj + geom_point(size=1.25) +  geom_line(size=0.25) + xlab("Target FDR level") 

				if (Mea == 'FDR') {
					obj <- obj +  ylab('False Discovery Proportion') +
						#	geom_hline(yintercept = 0, lty = 2, size = 0.25) +
							geom_abline(slope = 1, intercept=0,  size=0.25) 

				} else {
					obj <- obj + 
						#	facet_grid(SignalDensity ~ OrderPrior, scales='free') + 
							ylab('Number of Hits')
				}
				
				obj <- obj +
	
						scale_colour_manual(values=cols2) +
						scale_linetype_manual(values=ltys2) +
						scale_shape_manual(values=shapes2) +
						theme_bw() 

				plot.list[[Mea]] <- obj
				
			}


		}
		obj1 <- plot.list[['FDR']] + theme_bw(base_size = 14) 
		obj2 <- plot.list[['Power']]  + theme_bw(base_size = 14) 
		pdf(paste0(prefix, "_combined_curve.r2.pdf"), width=5, height=7)
		obj <- plot_grid(obj1, obj2, labels = "AUTO", ncol = 1, align = 'v', label_size=20)
		print(obj)
		dev.off()
	}
	
	
	cat('Finished!\n')
}
