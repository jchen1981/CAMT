require(qvalue)
require(IHW)
require(FDRreg)
require(adaptMT)
require(CAMT)
require(splines)
require(swfdr)
source('https://www.stat.uchicago.edu/~rina/sabha/All_q_est_functions.R')

# First set up working directory
source('Runtime/Code/run.time.func.R')
setwd('Runtime/Result')

cat('Start ...\n')
runtime.calc(methods = c('IHW', 'SABHA', 'BL', 'CAMT', 'FDRreg', 'AdaPT'), ns = c(1000, 4000, 16000, 64000, 256000, 1024000), iters = 1:3)
cat('Completed!\n')


