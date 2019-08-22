
clsapply <- function(dat, func, df=NULL, tempdir="~/project/tempZ",
		queque="1-hour", timespan="540", mem="8G", req.pack='stats',
		tidy=F) {
	# Apply to the row of dat by default (dim=2 for column)
	# Create temp dir
	if (!file.exists(tempdir)) {
		dir.create(tempdir)
	}
	setwd(tempdir)
	
	save(dat, func, df, file=("dat.RData"))
	
	nJob <- length(dat)
	# Create R script
	txt <- paste('
					args=(commandArgs(TRUE))
					if(length(args)==0){
						print("No arguments supplied.")
						##supply default values
					} else{
						for(i in 1:length(args)){
						eval(parse(text=args[[i]]))
					}
				}',
			paste(paste0('\nrequire(',req.pack, ')'), collapse="\n"),    
			'\n date()
					load("dat.RData")	
					if (is.list(dat)) {
						item <- dat[[part]]
					} 
					if (is.vector(dat)) {
						item <- dat[part]
					}           
					res0 <- func(item, df)
					save(res0, file=paste(part, ".res", sep=""))
					date()
			')
	writeLines(txt, "script.R")
	
	# Submit job
	cat("Submit jobs ...")
	prefix <- paste(c("J", sample(LETTERS, 4, repl=TRUE)), collapse="")
	for (part in 1:nJob) {
		rfile <- "script.R"
		rout <- paste(part, ".Rout", sep="")
			
		sh <- paste(
				"qsub",
				paste0("-N ", prefix, part),
				"-j y",
				"-cwd",
	#			"-t", timespan,
				"-q", queque,
				"-m abe",
				paste0("-l h_vmem=", mem),
				"-V",
				paste("-o ",  part, ".out", sep=""),
				paste("-e ",  part, ".err", sep=""),
				"-b y",
				paste("\'R CMD BATCH --no-save --no-restore \"--args  part=", part, "\" ", 
						rfile, " ", rout, "\'", sep="")
		)
		print(sh)
		system(sh)
	}
#	sbatch -J RP -t 240 -p general --mem 64000 -o RP.out -e RP.err --wrap='R CMD BATCH --no-save --no-restore ROC_SVA_PCA_ResultProcess_alpha.R RP.Rout'
	cat("Please wait ...\n")
	# To see if all the jobs have been finished
	while (TRUE) {
		Sys.sleep(20)
		output1 <- system("qstat ", intern=TRUE)
		output2 <- system("qstat ", intern=TRUE)
		missingfile <- NULL
		if (length(grep(prefix, output1)) == 0 & length(grep(prefix, output2)) == 0) {	
			cat("All jobs finished and begin to combine\n")
			res <- NULL
			for (part in 1:nJob) {
				resfile <- paste(part, ".res", sep="")
				if (file.exists(resfile)) {
					load(resfile)
					res[[paste(part)]] <-  res0
					rm(res0)
					cat(".")
				} else {
					# cat("\n", resfile, " not found!\n")
					missingfile <- c(missingfile, resfile)
				}
			}
			cat("\n")
			break
		}
		gc()
	}
	
	cat("Missing files:\n")
	if (length(missingfile) == 0) {
		cat("No jobs dropped!\n")
	} else{
		cat(paste(missingfile, collapse="\n"))
	} 
	# Clean
	if (tidy == TRUE) {
		cat("\nBegin to tidy up ...")
		system("rm *")
	}
	cat("\nDone!\n")
	return(res)
}


