require(ggplot2)
require(reshape)
require(cowplot)

# First set up working directory
###################################################################################################
# Plot - Annonation, theme
algs <- c('Oracle', 'CAMT',  'AdaPT',  'IHW', 'FDRreg', 'FDRregE', 'SABHA', 'BL', 'ST',  'BH')
cols <- c('Oracle'="#999999",'CAMT'="red",  'AdaPT'= 'blue','IHW'="#D55E00", 'FDRreg'="#CC79A7", 
		'FDRregE'="#CC79A7", 'SABHA'="#009E73", 'BL' = "#0072B2", 'BH'="gray60", 'ST'="gray30")
ltys <-  c('Oracle'=1, 'CAMT'=1, 'AdaPT'= 2,'IHW'= 1,  'FDRreg'=2, 'FDRregE'=2, 'SABHA'=1, 'BL'=2,'BH'=2, 'ST'=1)
shapes <- c('Oracle'=17,'CAMT'=19,   'AdaPT'=3,  'IHW'=2, 'FDRreg'=0, 'FDRregE'=0,'SABHA'=4, 'BL'=6, 'BH'=15, 'ST'=5)
###################################################################################################
alphalist <- seq(0.01, 0.2, len = 20)
algs2 <- c('CAMT',  'IHW', 'FDRreg', 'AdaPT', 'SABHA', 'BL', 'BH', 'ST')

IDs <- c('bottomly', 'pasilla', 'airway',   'yeast')
ID.ann <- c('Bottomly', 'Pasilla', 'Airway',  'Yeast protein')
names(ID.ann) <- IDs

plot.list <- NULL

for (ID in IDs) {
	load(file = file.path('RealData', 'Result', paste0('', ID, '.RData')))
	
	data <- res
	colnames(data) <- paste(alphalist)
	res <- melt(data)
	colnames(res) <- c("Method", 'Cutoff', 'value')
	
	res <- subset(res, Method %in% algs2  & Cutoff <= 0.10)
	res <- droplevels(res)
	res2 <- res
	res2$Cutoff <- factor(res2$Cutoff)
	res2 <- res2[as.numeric(res2$Cutoff) %in% (c(1, 2, 4, 6, 8, 10) ), ]
	res2$Cutoff <- as.numeric(as.character(res2$Cutoff))
	res2$Method <- factor(res2$Method, levels = algs2)
	res$Method <- factor(res$Method, levels = algs2)
	
	pdf(file.path('RealData', 'Result', paste0(ID, '.pdf')), height = 4, width = 5)
	obj <- ggplot(res, aes(x = Cutoff, y = value, group = Method,  colour = Method, linetype = Method)) +
			geom_line(size=0.25) +
			geom_point(data=res2, aes(x = Cutoff, y = value, shape = Method), size = 1) +
			ylab('Number of rejections') +
			xlab('Target FDR level') +
			xlim(c(0, 0.1)) +
			theme_bw() +
			ggtitle(ID.ann[ID]) +
			scale_linetype_manual(values = ltys[algs2]) +
			scale_colour_manual(values = cols[algs2]) +
			scale_shape_manual(values = shapes[algs2])
	print(obj)
	dev.off()	
	plot.list[[ID]] <- obj
}



pdf(file.path('RealData', 'Result', 'Combinedplot.AdaPT.Paper.pdf'), height = 6, width = 8)
plot_grid(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], labels = "AUTO", ncol = 2, align = 'v')
dev.off()

##############################################################################################################
alphalist <- seq(0.01, 0.2, len=20)
algs2 <- c('CAMT',  'IHW', 'FDRreg', 'AdaPT', 'SABHA', 'BL', 'BH', 'ST')

IDs <- c('ewas', 'mwas')
ID.ann <- c('CHD EWAS', 'AmericanGut Sex MWAS')
names(ID.ann) <- IDs

plot.list <- NULL
for (ID in IDs) {
	load(file = file.path('RealData', 'Result', paste0(ID, '.RData')))
	
	data <- res
	colnames(data) <- paste(alphalist)
	res <- melt(data)
	colnames(res) <- c("Method", 'Cutoff', 'value')
	
	res <- subset(res, Method %in% algs2 & Cutoff <= 0.20)
	res <- droplevels(res)
	res2 <- res
	res2$Cutoff <- factor(res2$Cutoff)
	res2 <- res2[as.numeric(res2$Cutoff) %in% (c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20) ), ]
	res2$Cutoff <- as.numeric(as.character(res2$Cutoff))
	res2$Method <- factor(res2$Method, levels = algs2)
	res$Method <- factor(res$Method, levels = algs2)
	
	pdf(file.path('RealData', 'Result', paste0(ID, '.pdf')), height = 4, width = 5)
	obj <- ggplot(res, aes(x = Cutoff, y = value, group = Method,  colour = Method, linetype = Method)) +
			geom_line(size=0.25) +
			geom_point(data=res2, aes(x = Cutoff, y = value, shape = Method), size = 1) +
			ylab('Number of rejections') +
			xlab('Target FDR level') +
			xlim(c(0, 0.2)) +
			theme_bw() +
			ggtitle(ID.ann[ID]) +
			scale_linetype_manual(values = ltys) +
			scale_colour_manual(values = cols) +
			scale_shape_manual(values = shapes)
	print(obj)
	dev.off()	
	plot.list[[ID]] <- obj
}


pdf(file.path('RealData', 'Result', 'Combinedplot.ewas.mwas.pdf'), height = 3, width = 8)
plot_grid(plot.list[[1]], plot.list[[2]], labels = "AUTO", ncol = 2, align = 'v')
dev.off()








