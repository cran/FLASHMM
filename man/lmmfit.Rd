% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmmfit.R
\name{lmmfit}
\alias{lmmfit}
\title{Fitting Linear Mixed-effects Models}
\usage{
lmmfit(
  Y,
  X,
  Z,
  d,
  theta0 = NULL,
  nBlocks = ceiling(nrow(Y) * ncol(Y) * 1e-08),
  method = c("REML", "ML"),
  max.iter = 50,
  epsilon = 1e-05,
  output.cov = TRUE,
  output.RE = FALSE
)
}
\arguments{
\item{Y}{A features-by-samples matrix of responses (genes-by-cells matrix of gene expressions for scRNA-seq).}

\item{X}{A design matrix for fixed effects, with rows corresponding to the columns of Y.}

\item{Z}{A design matrix for random effects, with rows corresponding to the columns of Y. Z = [Z1, ..., Zk], and Zi, i=1,...,k, is the design matrix for the i-th type random factor.}

\item{d}{A vector of (m1,...,mk), mi = ncol(Zi), number of columns in Zi. m1 + ... + mk = ncol(Z), number of columns in Z.}

\item{theta0}{A vector of initial values of the variance components, (s1, ...,sk, s_(k+1)), si = sigma_i^2, the variance component of the i-th type random effects. s_(k+1) = sigma^2, the variance component of model residual error.}

\item{nBlocks}{Number of blocks, used for blocking a big data to reduce the storage required in computing.}

\item{method}{The REML with Fisher scoring (FS) iterative algorithm, REML-FS.}

\item{max.iter}{The maximal number of iterations for the iterative algorithm.}

\item{epsilon}{Positive convergence tolerance. If the absolute value of the first partial derivative of log likelihood is less than epsilon, the iterations converge.}

\item{output.cov}{If TRUE, output the covariance matrices for the estimated coefficients, which are needed for testing contrasts.}

\item{output.RE}{If TRUE, output the best linear unbiased prediction (BLUP) of the random effects.}
}
\value{
A list containing the following components:

dlogL First partial derivatives of log-likelihoods for each feature (gene).

logLik Maximum log-likelihoods (ML) or log-restricted-likelihood (REML) for each feature (gene).

niter Nmbers of iterations for each feature (gene).

coef A matrix of estimated coefficients (fixed effects), each column corresponds to a feature (gene) and each row one covariate.

se A matrix of the standard errors of the estimated coefficients.

t A matrix of t-values for the fixed effects, equal to coef/se.

df Degrees of freedom.

p A matrix of two-sided p-values for the fixed effects.

cov A array of covariance matrices of the estimated coefficients (fixed effects).

theta A matrix of the estimated variance components, each column corresponds to a feature (gene) and each row one variance component. The last row is the variance component of the residual error.

se.theta Standard errors of the estimated theta.

RE A matrix of the best linear unbiased prediction (BLUP) of random effects.
}
\description{
lmmfit, a wrapper function of lmm, fits linear mixed-effects models (LMM) by sample-level data. The LMM parameters are estimated by restricted maximum likelihood (REML) with Fisher scoring (FS) gradient descent algorithm.
}
\examples{
#Generate data: X, Y, and Z.
set.seed(2024)

n <- 1e3
m <- 10
Y <- matrix(rnorm(m*n), m, n)
rownames(Y) <- paste0("Gene", 1:nrow(Y))

trt <- sample(c("A", "B"), n, replace = TRUE)
X <- model.matrix(~ 0 + trt)

q <- 20
sam <- rep(NA, n)
sam[trt == "A"] <- paste0("A", sample.int(q/2, sum(trt == "A"), replace = TRUE))
sam[trt == "B"] <- paste0("B", sample.int(q/2, sum(trt == "B"), replace = TRUE))
Z <- model.matrix(~ 0 + sam)
d <- ncol(Z)

#Fit LMM by the cell-level data
fit <- lmmfit(Y, X, Z, d = d)
str(fit)

#Hypothesis testing
lmmtest(fit)
lmmtest(fit, index = 2)
lmmtest(fit, contrast = cbind("B-A" = c(-1, 1)))

}
