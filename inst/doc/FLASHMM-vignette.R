## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = TRUE, results = "hide", message = FALSE---------------------------
##Install FLASHMM from CRAN.
# install.packages("FLASHMM")
##Install FLASHMM from Github.
devtools::install_github("https://github.com/Baderlab/FLASHMM")

##Load the package.
library(FLASHMM)

## ----dataset, message = FALSE-------------------------------------------------
set.seed(2412)
dat <- simuRNAseq(nGenes = 50, nCells = 1000, nsam = 25, ncls = 4, ntrt = 2, nDEgenes = 6)
names(dat)

#counts and meta data
counts <- dat$counts
metadata <- dat$metadata
head(metadata)
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

##(2) LMM fitting

##Option 1: fit LMM by lmmfit using cell-level data.
fit <- lmmfit(Y, X, Z, d = d)

##Option 2: fit LMM by lmm using summary-level data.
##- Computing summary statistics
n <- nrow(X)
XX <- t(X)%*%X; XY <- t(Y%*%X)
ZX <- t(Z)%*%X; ZY <- t(Y%*%Z); ZZ <- t(Z)%*%Z
Ynorm <- rowSums(Y*Y)
##- Fitting LMM
fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)

identical(fit, fitss)

##(3) Hypothesis tests
##Testing coefficients (fixed effects)
test <- lmmtest(fit)
#head(test)

##The testing t-value and p-values are also provided in the LMM fit.
range(test - cbind(t(fit$coef), t(fit$t), t(fit$p)))
#fit$coef[, 1:4]
#fit$t[, 1:4]
fit$p[, 1:4]
##

##Testing contrasts
##We can make comparisons using contrasts. For example, 
##the effects of treatment B vs A in all clusters can be tested 
##using the contrast constructed as follows:
ct <- numeric(ncol(X))
names(ct) <- colnames(X)
index <- grep("B", colnames(X))
ct[index] <- 1/length(index)
ct
test <- lmmtest(fit, contrast = ct)
head(test)


## ----include = TRUE, echo = TRUE, message = FALSE-----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("ExperimentHub", quietly = TRUE))
    BiocManager::install("ExperimentHub")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    BiocManager::install("SingleCellExperiment")

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
colnames(X) <- gsub("cls", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))

##Design matrix for random effects
Z <- model.matrix(~ 0 + as.factor(sam), data = dat$meta)
##Dimension of random effects
d <- ncol(Z)


## ----LMM fitting, echo = TRUE, message = FALSE, warning = FALSE---------------
  
##Option 1: fit LMM by cell-level data.
max.iter <- 100
epsilon <- 1e-5
fit <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)
names(fit)
##

##Option 2: fit LMM by summary-level data.
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
index <- grep(":", rownames(fit$coef))
contrast[index, ] <- 1/length(index)

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
##Coefficients, t-values, and p-values for the genes specific to a cell-type.
index <- grep(":", rownames(fit$coef))
ce <- fit$coef[index, cvg]
tv <- fit$t[index, cvg]
pv <- fit$p[index, cvg]

out <- data.frame(
	gene = rep(colnames(ce), nrow(ce)), 
	cluster = rep(rownames(ce), each = ncol(ce)),
	coef = c(t(ce)), t = c(t(tv)), p = c(t(pv)))

##Adjust p-values by FDR.
out$FDR <- p.adjust(out$p, method = "fdr")

##The top DE genes
##The DE genes with FDR < 0.05
out <- out[order(out$p), ]
rownames(out) <- NULL
out[out$FDR <= 0.05, ]

sessionInfo()


