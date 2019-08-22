
generate.example.data <- function (N = 100000, pi0 = 0.99) {
	
	truth <- rbinom(N, prob = 1 - pi0, size = 1)
	x0 <- gl(2, N / 2)
	x1 <- sample(x0)
	x2 <- sample(x0)
	z <- rnorm(N)
	
	ind1 <- x1 == 1 & x2 == 1 & truth == 1
	ind2 <- x1 == 2 & x2 == 1 & truth == 1
	ind3 <- x1 == 1 & x2 == 2 & truth == 1
	ind4 <- x1 == 2 & x2 == 2 & truth == 1
	
	z1 <- rnorm(N, 0,  1.5)
	z2 <- rnorm(N, 0,  2.0)
	z3 <- rnorm(N, 0,  2.5)
	z4 <- rnorm(N, 0,  3.0)
	z[ind1] <- z1[ind1]
	z[ind2] <- z2[ind2]
	z[ind3] <- z3[ind3]
	z[ind4] <- z4[ind4]
	
	pval <-  2 * pnorm(abs(z), lower = FALSE)
	pval1 <- 2 * pnorm(abs(z1), lower = FALSE)[ind1]
	pval2 <- 2 * pnorm(abs(z2), lower = FALSE)[ind2]
	pval3 <- 2 * pnorm(abs(z3), lower = FALSE)[ind3]
	pval4 <- 2 * pnorm(abs(z4), lower = FALSE)[ind4]
	
	return(list(pval = pval, x1 = x1, x2 = x2, truth = truth, z = z, 
					ind1 = ind1, ind2 = ind2, ind3 = ind3, ind4 = ind4,
					pval1 = pval1, pval2 = pval2, pval3 = pval3, pval4 = pval4))
}

plot.result <- function (camt.obj, data.obj, ann = '') {
	
	k1 <- unique(camt.obj$k[data.obj$ind1])
	k2 <- unique(camt.obj$k[data.obj$ind2])
	k3 <- unique(camt.obj$k[data.obj$ind3])
	k4 <- unique(camt.obj$k[data.obj$ind4])
	
	ecdf(data.obj$pval1) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k1, 1), add=TRUE, col='red')
	mtext(ann)
	
	ecdf(data.obj$pval2) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k2, 1), add=TRUE, col='red')
	
	ecdf(data.obj$pval3) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k3, 1), add=TRUE, col='red')
	
	ecdf(data.obj$pval4) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k4, 1), add=TRUE, col='red')
}

# All H1 hypotheses without any estimation uncertainty of pi0
generate.baseline.plot <- function (N = 10000) {
	
	z1 <- rnorm(N, 0,  1.5)
	z2 <- rnorm(N, 0,  2.0)
	z3 <- rnorm(N, 0,  2.5)
	z4 <- rnorm(N, 0,  3.0)
	
	pval1 <- 2 * pnorm(abs(z1), lower = FALSE)
	pval2 <- 2 * pnorm(abs(z2), lower = FALSE)
	pval3 <- 2 * pnorm(abs(z3), lower = FALSE)
	pval4 <- 2 * pnorm(abs(z4), lower = FALSE)
	
	
	target <- function(k, pval) {
		-sum(log(1 - k) - k * log(pval))
	}
	
	k1 <- optimize(target, c(0, 1), pval = pval1)$minimum
	k2 <- optimize(target, c(0, 1), pval = pval2)$minimum
	k3 <- optimize(target, c(0, 1), pval = pval3)$minimum
	k4 <- optimize(target, c(0, 1), pval = pval4)$minimum
	
	ecdf(pval1) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k1, 1), add=TRUE, col='red')
	
	mtext('Baseline')
	
	ecdf(pval2) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k2, 2), add=TRUE, col='red')
	
	ecdf(pval3) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k3, 1), add=TRUE, col='red')
	
	ecdf(pval4) %>% plot(log='x', xlim = c(1e-12, 1), xlab = 'p value', main = '')
	curve(pbeta(x, 1 - k4, 1), add=TRUE, col='red')
	
}

require(CAMT)

# Set working directory
set.seed(123)
pdf('Illustration_f1.pdf', width = 12, height = 12)
par(mfrow = c(4, 4))


data.obj <- generate.example.data(N = 50000, pi0 = 0.99)
camt.obj <- camt.fdr(pvals = data.obj$pval, pi0.var = NULL, f1.var = cbind(data.obj$x1, data.obj$x2), pvals.cutoff = min(data.obj$pval))
plot.result(camt.obj, data.obj, 'pi0 = 0.99')

data.obj <- generate.example.data(N = 50000, pi0 = 0.95)
camt.obj <- camt.fdr(pvals = data.obj$pval, pi0.var = NULL, f1.var = cbind(data.obj$x1, data.obj$x2), pvals.cutoff = min(data.obj$pval))
plot.result(camt.obj, data.obj, 'pi0 = 0.95')

data.obj <- generate.example.data(N = 50000, pi0 = 0.80)
camt.obj <- camt.fdr(pvals = data.obj$pval, pi0.var = NULL, f1.var = cbind(data.obj$x1, data.obj$x2), pvals.cutoff = min(data.obj$pval))
plot.result(camt.obj, data.obj, 'pi0 = 0.80')

# For comparsion
generate.baseline.plot (N = 10000)

dev.off()



