## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = TRUE, results = "hide", message = FALSE---------------------------
##Install FLASHMM from CRAN.
# install.packages("FLASHMM")
##Install FLASHMM from Github.
# devtools::install_github("https://github.com/Baderlab/FLASHMM")

##Load the package.
library(FLASHMM)

## ----Reference dataset, echo = TRUE, message = FALSE--------------------------
##Generate a reference dataset by simuRNAseq function.
set.seed(2502)
refdata <- simuRNAseq(nGenes = 1000, nCells = 10000)

##counts
counts <- refdata$counts
##metadata
metadata <- refdata$metadata
head(metadata)

rm(refdata)

## ----Simulate dataset, echo = TRUE, message = FALSE---------------------------
##Generate the scRNA-seq data by simuRNAseq function.
set.seed(2503)
dat <- simuRNAseq(counts, metadata = metadata, nGenes = 100, nCells = 100000, nsam = 25, 
       ncls = 10, ntrt = 2, nDEgenes = 10, minbeta = 0.5, maxbeta = 1, var.randomeffects = 0.1)
str(dat)
rm(counts, metadata)

## ----Expression matrix, echo = TRUE, message = FALSE--------------------------
##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + dat$counts)
dat$counts <- NULL #Remove the counts.

##Design matrix for fixed effects
X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = dat$meta)

##Design matrix for random effects
Z <- model.matrix(~ 0 + as.factor(sam), data = dat$meta)
##Dimension of random effects
d <- ncol(Z)

rm(dat)

## ----LMM fitting, echo = TRUE, message = FALSE, warning = FALSE---------------
##Option 1: fit LMM by cell-level data.
max.iter <- 100
epsilon <- 1e-5
fit <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)

##Option 2: fit LMM by summary-level data.
##(1) Compute the summary-level data.
n <- nrow(X)
XX <- t(X)%*%X
XY <- t(Y%*%X)
ZX <- t(Z)%*%X
ZY <- t(Y%*%Z)
ZZ <- t(Z)%*%Z
Ynorm <- rowSums(Y*Y)

rm(X, Y, Z) #release the memory.

##(2) Fit LMM.
fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, 
             max.iter = max.iter, epsilon = epsilon)

identical(fit, fitss) #The two LMM fits are identical.
rm(fitss)


## -----------------------------------------------------------------------------
##Testing coefficients (fixed effects)
test <- lmmtest(fit)
#head(test)

##The testing t-value and p-values are also provided in the LMM fit.
range(test - cbind(t(fit$coef), t(fit$t), t(fit$p))) #identical
#fit$t[, 1:4]
#fit$p[, 1:4]
fit$coef[, 1:4]

## -----------------------------------------------------------------------------
contrast <- cbind("BvsA" = numeric(nrow(fit$coef)))
index <- grep(":", rownames(fit$coef))
contrast[index, ] <- 1/length(index)

##Test the contrast.
test <- lmmtest(fit, contrast = contrast)
head(test)


## -----------------------------------------------------------------------------
BvsA <- paste0(paste0("cls", 1:10, ":trtB"), collapse = "+")
BvsA <- paste0("(", BvsA, ")/10")
BvsA
contrast <- contrast.matrix(contrast = c(BvsA = BvsA), model.matrix.names = rownames(fit$coef))
test <- lmmtest(fit, contrast = contrast)
head(test)

## ----DE genes, echo = TRUE, message = FALSE, warning = FALSE------------------
##(1) Check the genes for which LMM fitting converges.
gc <- (apply(abs(fit$dlogL), 2, max) < epsilon) 
sum(gc) 

##(2) Hypothesis testing for variance components:
##    H0, theta </= 0 vs H1, theta > 0.
z <- fit$theta["var1", ]/fit$se.theta["var1", ]
p <- pnorm(z, lower.tail = FALSE)
sum(p < 0.05)

##(3) The DE genes specific to a cell-type
##Coefficients, t-values, and p-values for the genes specific to a cell-type.
index <- grep(":", rownames(fit$coef))
ce <- fit$coef[index, gc]
tv <- fit$t[index, gc]
pv <- fit$p[index, gc]

out <- data.frame(
	gene = rep(colnames(ce), nrow(ce)), 
	cluster = rep(rownames(ce), each = ncol(ce)),
	coef = c(t(ce)), t = c(t(tv)), p = c(t(pv)))

##Adjust p-values by FDR.
out$FDR <- p.adjust(out$p, method = "fdr")

##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
out <- out[order(out$p), ]
rownames(out) <- NULL
out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]

##The true DE genes
#dat$DEgenes

## ----LMM_ML, echo = TRUE, message = FALSE, warning = FALSE--------------------
##Fitting LMM by ML method
#fit <- lmmfit(Y, X, Z, d = d, method = "ML", max.iter = max.iter, epsilon = epsilon)
fit <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, method = "ML",
             max.iter = max.iter, epsilon = epsilon)

##The DE genes specific to a cell-type
##Coefficients, t-values, and p-values
index <- NULL
index <- grep(":", rownames(fit$coef))
ce <- fit$coef[index, gc]
tv <- fit$t[index, gc]
pv <- fit$p[index, gc]

out <- NULL
out <- data.frame(
	gene = rep(colnames(ce), nrow(ce)), 
	cluster = rep(rownames(ce), each = ncol(ce)),
	coef = c(t(ce)), t = c(t(tv)), p = c(t(pv)))

##Adjusting p-values by FDR
out$FDR <- p.adjust(out$p, method = "fdr")

##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
out <- out[order(out$p), ]
rownames(out) <- NULL
out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]

##
sessionInfo()


