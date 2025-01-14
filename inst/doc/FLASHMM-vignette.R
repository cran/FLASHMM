## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = TRUE, results = "hide", message = FALSE---------------------------
#Install FLASHMM from Github.
devtools::install_github("https://github.com/Baderlab/FLASHMM")

#Load the package.
library(FLASHMM)

## ----dataset, message = FALSE-------------------------------------------------
set.seed(2412)
dat <- simuRNAseq(nGenes = 50, nCells = 1000, 
                  nsam = 25, ncls = 4, ntrt = 2, nDEgenes = 6)

names(dat)
##

#counts and meta data
counts <- dat$counts
metadata <- dat$metadata
rm(dat)

## -----------------------------------------------------------------------------
##(1) Model design
##Y: gene expression profile (log-transformed counts)
##X: design matrix for fixed effects
##Z: design matrix for random effects

Y <- log(counts + 1) 
X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = metadata)
Z <- model.matrix(~ 0 + sam, data = metadata)
d <- ncol(Z)
##

##(2) LMM fitting
##Fit LMM by lmmfit using cell-level data.
fit <- lmmfit(Y, X, Z, d = d)

##Fit LMM by lmm using summary-level data computed as follows.
##- Computing summary statistics
n <- nrow(X)
XX <- t(X)%*%X; XY <- t(Y%*%X)
ZX <- t(Z)%*%X; ZY <- t(Y%*%Z); ZZ <- t(Z)%*%Z
Ynorm <- rowSums(Y*Y)
##- Fitting LMM
fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)

identical(fit, fitss)
##

##Fit LMM by lmm using summary-level data computed by sslmm.
##- Computing summary statistics
ss <- sslmm(X, Y, Z)
##- Fitting LMM
fitss <- lmm(summary.stats = ss, d = d)

identical(fit, fitss)
##

##(3) Hypothesis tests
test <- lmmtest(fit)
#head(test)

##t-values
all(t(fit$t) == test[, grep("_t", colnames(test))])
fit$t[, 1:5]
##

##p-values
all(t(fit$p) == test[, grep("_p", colnames(test))])
fit$p[, 1:5]

## ----Reference dataset, echo = TRUE, message = FALSE--------------------------
library(ExperimentHub)

##Load data.
eh <- ExperimentHub()
#query(eh, "Kang")
sce <- eh[["EH2259"]]

##Remove undetected genes.
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
##Remove cells with few or many detected genes (outliers).
nGenes <- colSums(counts(sce) > 0)
bx <- boxplot(log(nGenes), plot = FALSE)
sce <- sce[, nGenes >= exp(bx$stats[1]) & nGenes <= exp(bx$stats[5])]
##Remove lowly expressed genes.
##counts per cell (cpc) 
cpc <- rowSums(counts(sce))/ncol(sce)
sce <- sce[(rowSums(counts(sce) > 1) >= 10) & (cpc > 0.005), ]

##counts and metadata
counts <- assay(sce, "counts")
coldata <- as.data.frame(colData(sce))
head(coldata)
all(colnames(counts) == rownames(coldata))
dim(counts)
rm(eh, sce)


## ----Simulate dataset, echo = TRUE, message = FALSE---------------------------
##Specify which columns represent samples, treatments, and cell-types.
colnames(coldata)[colnames(coldata) == "ind"] <- "sam"  #samples
colnames(coldata)[colnames(coldata) == "stim"] <- "trt" #treatments
colnames(coldata)[colnames(coldata) == "cell"] <- "cls" #cell-types
coldata <- coldata[, c("sam", "trt", "cls")]
head(coldata)
##

##Generate the dataset by simuRNAseq function.
set.seed(2412)
dat <- simuRNAseq(counts, nGenes = 100, nCells = 120000, metadata = coldata, 
                  nsam = 25, ncls = 10, ntrt = 2, nDEgenes = 10, 
                  minbeta = 0.5, maxbeta = 1, var.randomeffects = 0.1)
str(dat)

##Remove the reference dataset that is no longer needed in the following analysis.
rm(counts, coldata)

## ----Expression matrix, echo = TRUE, message = FALSE--------------------------
##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + dat$counts)
dat$counts <- NULL #Remove the counts.

##Design matrix for fixed effects
X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = dat$meta)
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))

##Design matrix for random effects
Z <- model.matrix(~ 0 + as.factor(sam), data = dat$meta)
##Dimension of random effects
d <- ncol(Z)


## ----LMM fitting, echo = TRUE, message = FALSE, warning = FALSE---------------
  
##(1) Fit LMM by cell-level data.
max.iter <- 100
epsilon <- 1e-5
fit <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)
str(fit)
##

##(2) Fit LMM by summary-level data.
##- Compute the summary-level data.
n <- nrow(X)
XX <- t(X)%*%X
XY <- t(Y%*%X)
ZX <- t(Z)%*%X
ZY <- t(Y%*%Z)
ZZ <- t(Z)%*%Z
Ynorm <- rowSums(Y*Y)

rm(X, Y, Z) #release the memory.

##- Fit LMM.
fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, 
             max.iter = max.iter, epsilon = epsilon)
identical(fit, fitss)
rm(fitss)

##Test the treatment effect within all cell-types.
##- Construct a contrast by summing the treatment effects across cell-types.
contrast <- cbind("trt" = numeric(nrow(fit$coef)))
contrast[grep(":", rownames(fit$coef)), ] <- 1

##- Test the contrast.
test <- lmmtest(fit, contrast = contrast)
head(test)


## ----DE genes, echo = TRUE, message = FALSE, warning = FALSE------------------

##(1) Check which LMM fittings converge.
cvg <- (apply(abs(fit$dlogL), 2, max) < epsilon) 
sum(cvg) 
##

##(2) Hypothesis testing for variance components:
##    H0, theta </= 0 vs H1, theta > 0.
z <- fit$theta["var1", ]/fit$se.theta["var1", ]
p <- pnorm(z, lower.tail = FALSE)
sum(p <= 0.05)
##

##(3) The DE genes specific to a cell-type
##Coefficients and p-values for the genes specific to a cell-type.
index <- grep(":", rownames(fit$coef))
beta <- t(fit$coef[index, cvg])
p <- t(fit$p[index, cvg])

##Adjust p-values by FDR.
padj <- matrix(p.adjust(c(p), method = "fdr"), nrow = nrow(p), ncol = ncol(p))

##The DE genes specific to a cell cluster with FDR < 0.05.
degenes <- NULL
for (j in 1:ncol(p)){
  i <- (padj[, j] < 0.05)
  if (sum(i) > 0) degenes <- rbind(degenes, data.frame(gene = rownames(p)[i], cluster = j, coef = beta[i, j], p = p[i, j], FDR = padj[i, j]))
}
rownames(degenes) <- NULL
degenes

##The simulated DE genes
dat$DEgenes

sessionInfo()


