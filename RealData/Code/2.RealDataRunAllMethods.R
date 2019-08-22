require(qvalue)
require(IHW)
require(FDRreg)
require(adaptMT)
require(CAMT)
require(splines)
require(swfdr)
source('https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R')

# First set up working directory
source('RealData/Code/run.real.func.R')

for (i in paste(1:6)) {
	
	study <- switch(i, '1' = 'bottomly', '2' = 'pasilla', '3' = 'airway', '4' = 'yeast', '5' = 'ewas', '6' = 'mwas')
	
	cat(study, "\n")
	cat("Load data ...\n")
	
	data <- readRDS(file.path('RealData', 'Data', paste0(study, '.p.value.rds')))
	res <- run.all.method(data[, "pvalue"], data[, "covariate"])
	save(res, file=file.path('RealData', 'Result', paste(study, ".Rdata", sep = "")))
	cat("Completed!\n")
}