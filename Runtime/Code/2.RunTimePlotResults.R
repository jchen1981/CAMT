# First set up working directory
setwd('Runtime/Result')

#####################################################################
# Load data
methods <- c('CAMT', 'AdaPT', 'IHW', 'FDRreg', 'SABHA', 'BL')

ns <- 1000 * 4 ^ (0:5)
iters <- 1:3
res <- array(NA, c(length(methods), length(ns), length(iters)), 
		dimnames=list(Method = methods, N = paste(ns), I = paste(iters)))
for (n in ns) {
	for (method in methods) {
		for (i in iters) {
			file.name <- paste0(method, '.', n, '.', i, '.txt')
			if(file.exists(file.name)) {
				res[method, paste(n), paste(i)] <- scan(file=file.name)
			}
			
		}
	}
}

#####################################################################
# Plot data
require(ggplot2)
require(reshape2)
require(ggthemes)
library(scales) 

data <- melt(res)
col.algs2  <- c('CAMT'="red", 'AdaPT'='blue','IHW'="#D55E00", 'FDRreg'="#CC79A7", 'SABHA'="#009E73", 'BL'="#0072B2")
line.algs2 <- c('CAMT'=1, 'AdaPT'= 2,'IHW'= 1,  'FDRreg'=2, 'SABHA'=1, 'BL'=2)
shape.algs2 <- c('CAMT'=19,  'AdaPT'=3, 'IHW'=2, 'FDRreg'=0, 'SABHA'=4, 'BL'=6)


pdf('RunTimeComparsion.pdf', width = 5, height = 4)
data2 <- aggregate(value ~ Method + N, data, mean)
data2 <- subset(data2, Method %in% c('CAMT',  'IHW', 'FDRreg', 'AdaPT', 'SABHA', 'BL'))
data2$Method <- factor(data2$Method, levels = c('CAMT',  'IHW', 'FDRreg', 'AdaPT', 'SABHA', 'BL'))
levels(data2$Method) <- c('CAMT', 'IHW', 'FDRreg', 'AdaPT', 'SABHA', 'BL')
obj <- ggplot(data2, aes(x=N, y=value, colour=Method, linetype=Method, shape=Method)) +
		geom_point(size=2) +
		geom_line(size=0.25) +
		ylab('Seconds') +
		xlab('Sample size') +
		scale_colour_manual(values=col.algs2) +
		scale_linetype_manual(values=line.algs2) +
		scale_shape_manual(values=shape.algs2) +
		scale_x_continuous(trans = log10_trans(),
				breaks = trans_breaks("log10", function(x) 10^x),
				labels = trans_format("log10", math_format(10^.x))) +
		theme_bw()
print(obj)
dev.off()

