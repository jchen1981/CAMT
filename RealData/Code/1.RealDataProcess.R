# First set up the working directory
setwd('./RealData/Data')

######################################################################################################################
# Bottomly
# Bottomly, Daniel, et al. "Evaluating gene expression in C57BL/6J and DBA/2J mouse striatum using RNA-Seq and 
# microarrays." PloS one 6.3 (2011): e17820
######################################################################################################################
library("IHWpaper")
bottomly <- analyze_dataset("bottomly")

pvals <- bottomly$pvalue
ind <- which(!is.na(pvals))
pvals <- pvals[ind]
x1 <- log(bottomly$baseMean)
x1 <- x1[ind]

reorder <- rev(order(x1))
x1 <- x1[reorder]
pvals <- pvals[reorder]

hist(pvals)
length(pvals)

bottomly <- cbind(pvalue = pvals, covariate = x1)
saveRDS(bottomly, file='bottomly.p.value.rds')

######################################################################################################################
# Pasilla
# Brooks, Angela N., et al. "Conservation of an RNA regulatory map between Drosophila and mammals." Genome research 
# 21.2 (2011): 193-202.
######################################################################################################################
library("DESeq")
library("pasilla")
data("pasillaGenes")

cds <- estimateSizeFactors(pasillaGenes)
cds <- estimateDispersions(cds)
fit1 <- fitNbinomGLMs(cds, count ~ type + condition)
fit0 <- fitNbinomGLMs(cds, count ~ type)
res <- data.frame(
		filterstat = rowMeans(counts(cds)),
		pvalue = nbinomGLMTest(fit1, fit0),
		row.names = featureNames(cds))
ind <- which(!is.na(res$pvalue))
res <- res[ind, ]

pvals <- res$pvalue
x1 <- log(res[, 1])

reorder <- rev(order(x1))
x1 <- x1[reorder]
pvals <- pvals[reorder]

hist(pvals)
length(pvals)

pasilla <- cbind(pvalue = pvals, covariate = x1)
saveRDS(pasilla, 'pasilla.p.value.rds')

######################################################################################################################
# Airway
# Himes, Blanca E., et al. "RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene
# that modulates cytokine function in airway smooth muscle cells." PloS one 9.6 (2014): e99625.
######################################################################################################################
library("DESeq2")
library("dplyr")

data("airway", package = "airway")
dds <- DESeqDataSet(se = airway, design = ~ cell + dex) %>% DESeq
deRes <- as.data.frame(results(dds))

pvals <- deRes$pvalue
ind <- which(!is.na(pvals))
pvals <- pvals[ind]
x <- log(deRes$baseMean)
x <- x[ind]

reorder <- rev(order(x))
x1 <- x[reorder]
pvals <- pvals[reorder]

hist(pvals)
length(pvals)
airway <- cbind(pvalue = pvals, covariate = x1)

saveRDS(airway, 'airway.p.value.rds')

######################################################################################################################
# Yeast
# Dephoure, Noah, and Steven P. Gygi. "Hyperplexing: a method for higher-order multiplexed quantitative proteomics 
# provides a map of the dynamic response to rapamycin in yeast." Sci. Signal. 5.217 (2012): rs2-rs2.
######################################################################################################################
library("DESeq2")
library("dplyr")
library("IHWpaper")

proteomics_file <- system.file(
		"extdata/real_data",
		"science_signaling.csv",
		package = "IHWpaper"
)

proteomics_df <- read.csv(proteomics_file, stringsAsFactors = F)

proteomics_df$pvalue <- rank(
		proteomics_df$p1,
		ties.method="first"
) * proteomics_df$p1 / nrow(proteomics_df) 

pvals <- proteomics_df$pvalue
x <- log(proteomics_df$X..peptides)

reorder <- rev(order(x))
x1 <- x[reorder]
pvals <- pvals[reorder]

hist(pvals)
length(pvals)

yeast <- cbind(pvalue = pvals, covariate = x1)
saveRDS(yeast, 'yeast.p.value.rds')


######################################################################################################################
# EWAS - data directly from the paper, included in the package
# Wijnands, Kim PJ, et al. "Genome-wide methylation analysis identifies novel CpG loci for perimembranous ventricular
# septal defects in human." Epigenomics 9.3 (2017): 241-251.
######################################################################################################################
#pvals <- Haven.df[, 'pvalue']
#x1 <- Haven.df[, 'mean']
#
#reorder <- rev(order(x1))
#x1 <- x1[reorder]
#pvals <- pvals[reorder]
#
#ewas <- cbind(pvalue = pvals, covariate = x1)
#rownames(ewas) <- rownames(Haven.df)[reorder]
#
#saveRDS(ewas, 'ewas.p.value.rds')

######################################################################################################################
# MWAS - data retrieved from figshare biom file "deblur_125nt_no_blooms_normed.biom" - doi:10.6084/m9.figshare.6137198
# McDonald, Daniel, et al. "American gut: an open platform for citizen science microbiome research." mSystems (2018) :
# e00031-18.
######################################################################################################################
# The downloaded data were already normalized
# Rare OTUs occured in less than or equal to 20 subjects were excluded (>0.2% prevalence) in the data
meta.dat <- readRDS(file='MWAS/amgut.meta.dat.rds')
otu.tab <- readRDS(file='MWAS/amgut.otu.dat.rds')
otu.name <- readRDS(file='MWAS/amgut.otu.name.rds')

# We select subjects from United States and adults 
ind <- meta.dat$bmi_cat %in% c('Normal') & meta.dat$country_residence %in% c('United States') & 
		meta.dat$age_cat %in% c('teen', '20s', '30s', '40s', '50s', '60s') & meta.dat$sex %in% c('female', 'male')

otu.tab <- otu.tab[, ind]
meta.dat <- meta.dat[ind, ]
sex <- meta.dat$sex
sex <- factor(sex)
meta.dat <- droplevels(sex)

# Further discard OTU occuring in less than 5 subjects 
ind <- rowSums(otu.tab != 0) >= 5
otu.tab <- otu.tab[ind, ]
otu.name <- otu.name[ind, ]

# Wilcox rank sum test
pvals <- apply(otu.tab, 1, function (x) wilcox.test(x ~ sex)$p.value)
x <- rowSums(otu.tab != 0)

reorder <- rev(order(x))
x1 <- x[reorder]
pvals <- pvals[reorder]

hist(pvals)
length(pvals)

mwas <- data.frame(pvalue = pvals, covariate = x1, otu.name[reorder, ])
rownames(mwas) <- rownames(otu.tab)[reorder]

saveRDS(mwas, 'mwas.p.value.rds')

# Compare to traditional filtering based method
otu.tab2 <- otu.tab[rowMeans(otu.tab != 0) >= 0.1, ]
pvals2 <- apply(otu.tab2, 1, function (x) wilcox.test(x ~ sex)$p.value)
sum(qvalue(pvals2)$qvalue <= 0.1)  # 116
sum(p.adjust(pvals2, 'fdr') <= 0.1) # 85

otu.tab2 <- otu.tab[rowMeans(otu.tab != 0) >= 0.2, ]
pvals2 <- apply(otu.tab2, 1, function (x) wilcox.test(x ~ sex)$p.value)
sum(qvalue(pvals2)$qvalue <= 0.1)  # 71
sum(p.adjust(pvals2, 'fdr') <= 0.1) # 50
######################################################################################################################