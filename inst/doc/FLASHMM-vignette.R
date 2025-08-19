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
refdata <- simuRNAseq(nGenes = 100, nCells = 10000)
counts <- refdata$counts #counts
metadata <- refdata$metadata #metadata
head(metadata)

rm(refdata)

## ----Simulate dataset, echo = TRUE, message = FALSE---------------------------
##Generate the scRNA-seq data by simuRNAseq function.
set.seed(2503)
dat <- simuRNAseq(counts, metadata = metadata, nGenes = 100, nCells = 100000, 
       nsam = 25, ncls = 10, ntrt = 2, nDEgenes = 10)
str(dat)

##Samples (subjects) nested in one treatment A or B
#table(dat$metadata$sam)
#table(dat$metadata$trt)
all(grep("A", dat$metadata$sam) %in% which(dat$metadata$trt == "A"))
all(grep("B", dat$metadata$sam) %in% which(dat$metadata$trt == "B"))

rm(counts, metadata) #releasing memory

## ----Expression matrix, echo = TRUE, message = FALSE, tidy = TRUE, tidy.opts = list(width.cutoff = 90)----

##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + dat$counts)

dat$counts <- NULL #releasing memory

##Since the simulated data contains no covariates other than libsize (library size) and cls (clsters or cell types), the design matrix for fixed effects is given as follows. 
X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = dat$metadata)

##Note that the samples in the simulated data are nested in one treatment A or B. Design matrix for single-component random effects is given by
Z <- model.matrix(~ 0 + as.factor(sam), data = dat$metadata)
d <- ncol(Z)

##Suppose the data contains different measurement time points within a sample, denoted as 'time', which are randomly generated. To account for the variability of the time within a sample, we may use the two-component random effects design matrix: 
set.seed(250801)
n <- nrow(X)
dat$metadata$time <- runif(n, 1, 1.5)*sample(1:2, n, replace = TRUE)
Za <- model.matrix(~ 0 + sam + sam:time, data = dat$metadata)
da <- c(ncol(Za)/2, ncol(Za)/2) #dimension
range(Za[, 1:d] - Z) #identical to the single-component design

rm(dat) #releasing memory

## ----LMM fitting, echo = TRUE, message = FALSE, warning = FALSE---------------
##Option 1: Fit LMM based on cell-level data.
max.iter <- 100
epsilon <- 1e-5

##Fit the LMM with one-component random effects.
fit1 <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)

##Fit the LMM with Two-component random effects.
fit2 <- lmmfit(Y, X, Za, d = da, max.iter = max.iter, epsilon = epsilon)

##Option 2: Fit LMM based on summary-level data.
##(1) Compute the summary-level data.
n <- nrow(X)
XX <- t(X)%*%X
XY <- t(Y%*%X)
ZX <- t(Z)%*%X
ZY <- t(Y%*%Z)
ZZ <- t(Z)%*%Z
Ynorm <- rowSums(Y*Y)

##The summary-level data can also be computed by sslmm function:
ss <- sslmm(X, Y, Z)

XX <- ss$XX; XY <- ss$XY
ZX <- ss$ZX; ZY <- ss$ZY; ZZ <- ss$ZZ
n <- ss$n; Ynorm <- ss$Ynorm

rm(X, Y, Z, Za) #releasing memory.

##(2) Fit LMM using lmm function.
fitss <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, 
             max.iter = max.iter, epsilon = epsilon)

##The two LMM fits are identical.
identical(fit1, fitss) 
rm(fitss)
str(fit1)

## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
##Check the genes for which LMM fitting converges.
convg <- (apply(abs(fit1$dlogL), 2, max) < epsilon) 
sum(convg) 

## -----------------------------------------------------------------------------
##Testing coefficients (fixed effects)
test <- lmmtest(fit1)
##The t-value and p-values are identical with those provided in the LMM fit.
range(test - cbind(t(fit1$coef), t(fit1$t), t(fit1$p)))

fit1$coef[, 1:4]

#fit1$t[, 1:4]
#fit1$p[, 1:4]

## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
##The DE genes specific to a cell-type
##Coefficients, t-values, and p-values for the genes specific to a cell-type.
index <- grep(":", rownames(fit1$coef))
ce <- fit1$coef[index, ]
tv <- fit1$t[index, ]
pv <- fit1$p[index, ]

out <- data.frame(
	gene = rep(colnames(ce), nrow(ce)), 
	cluster = rep(rownames(ce), each = ncol(ce)),
	coef = c(t(ce)), t = c(t(tv)), p = c(t(pv)))

##FDR.
out$FDR <- p.adjust(out$p, method = "fdr")

##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
out <- out[order(out$p), ]
rownames(out) <- NULL
out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]


## -----------------------------------------------------------------------------
##Construct the contrast.
contrast <- cbind("BvsA" = numeric(nrow(fit1$coef)))
index <- grep(":", rownames(fit1$coef))
contrast[index, ] <- 1/length(index)
length(index)

##Test the contrast.
test <- lmmtest(fit1, contrast = contrast)
head(test)

## -----------------------------------------------------------------------------
ncls <- 10 #length(index)
sumeff <- paste0(paste0("cls", 1:ncls, ":trtB"), collapse = "+")
sumeff <- paste0("(", sumeff, ")/", ncls)
#sumeff
contrast <- contrast.matrix(
            contrast = c(BvsA = sumeff), 
            model.matrix.names = rownames(fit1$coef))
test <- lmmtest(fit1, contrast = contrast)
head(test)

## -----------------------------------------------------------------------------
##Using z-test for testing variance components
i <- grep("var1", rownames(fit1$theta))  #The (first) variance component
z <- fit1$theta[i, ]/fit1$se.theta[i, ]   #z-statistics

##One-sided z-test p-values for hypotheses:
##H0: theta <= 0 vs H1: theta > 0
p <- pnorm(z, lower.tail = FALSE)

##Number of significant p-values
sum(p < 0.05/length(p)) #Bonferroni correction


## -----------------------------------------------------------------------------
##(1) z-test for testing the second variance component
##Z-statistics for the second variance component
i <- grep("var2", rownames(fit2$theta)) 
z <- fit2$theta[i, ]/fit2$se.theta[i, ] 
##One-sided z-test p-values for hypotheses:
##H0: theta <= 0 vs H1: theta > 0
p <- pnorm(z, lower.tail = FALSE)

##number of significant p-values
sum(p.adjust(p, method = "fdr") < 0.05)
sum(p < 0.05/length(p)) #Bonferroni correction

##(2) LRT for testing the second variance component
LRT <- 2*(fit2$logLik - fit1$logLik)
pLRT <- pchisq(LRT, df = 1, lower.tail = FALSE)

##number of significant p-values
sum(p.adjust(pLRT, method = "fdr") < 0.05)
sum(pLRT < 0.05/length(pLRT)) #Bonferroni correction

##QQ-plot
qqplot(runif(length(p)), p, xlab = "Uniform quantile", ylab = "Z-test p-value")
abline(0, 1, col = "gray")
qqplot(runif(length(pLRT)), pLRT, xlab = "Uniform quantile", ylab = "LRT p-value")
abline(0, 1, col = "gray")

## ----LMM_ML, echo = TRUE, message = FALSE, warning = FALSE--------------------
##Fitting LMM using ML method
#fit3 <- lmmfit(Y, X, Z, d = d, method = "ML", max.iter = max.iter, epsilon = epsilon)
fit3 <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, method = "ML",
             max.iter = max.iter, epsilon = epsilon)

##The DE genes specific to a cell-type
##Coefficients, t-values, and p-values
index <- grep(":", rownames(fit3$coef))
ce <- fit3$coef[index, ]
tv <- fit3$t[index, ]
pv <- fit3$p[index, ]
out <- data.frame(
	gene = rep(colnames(ce), nrow(ce)), 
	cluster = rep(rownames(ce), each = ncol(ce)),
	coef = c(t(ce)), t = c(t(tv)), p = c(t(pv)))

##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
out$FDR <- p.adjust(out$p, method = "fdr") #FDR
out <- out[order(out$p), ]
rownames(out) <- NULL
out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]

## ----echo = TRUE, message = TRUE----------------------------------------------
sessionInfo()

