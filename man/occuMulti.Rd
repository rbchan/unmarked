\name{occuMulti}

\alias{occuMulti}

\title{Fit the Rota et al. (2016) Multi-species Occupancy Model}

\usage{occuMulti(detformulas, stateformulas, data, maxOrder, penalty=0, boot=30, 
       starts, method="BFGS", se=TRUE, engine=c("C","R"), silent=FALSE, ...)}

\arguments{
    \item{detformulas}{Character vector of formulas for the detection models, one per species.}
    \item{stateformulas}{Character vector of formulas for the natural parameters. 
      To fix a natural parameter at 0, specify the corresponding formula as \code{"0"} or \code{"~0"}.}
    \item{data}{An \code{\link{unmarkedFrameOccuMulti}} object}
    \item{maxOrder}{Optional; specify maximum interaction order. Defaults to 
      number of species (all possible interactions). Reducing this value may
      speed up optimization if you aren't interested in higher-order interactions.}
    \item{penalty}{Penalty term for likelihood. The total penalty is calculated
      as \code{penalty * 0.5 * sum(paramvals^2)}. Defaults to 0 (no penalty).}
    \item{boot}{Number of bootstrap samples to use to generate the 
      variance-covariance matrix when \code{penalty > 0}.}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{Logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Either "C" to use fast C++ code or "R" to use native R
      code during the optimization.}
    \item{silent}{Boolean; if \code{TRUE}, suppress warnings.}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{This function fits the multispecies occupancy model of Rota et al (2016).}

\details{

See \code{\link{unmarkedFrameOccuMulti}} for a description of how to supply data 
to the \code{data} argument.

\code{occuMulti} fits the multispecies occupancy model from Rota et al. (2016),
for two or more interacting species.
The model generalizes the standard single-species occupancy model from
MacKenzie et al. (2002).
The latent occupancy state at site \eqn{i} for a set of \eqn{s} potentially
interacting species is a vector \eqn{\mathbf{Z}_i}{Z_i} of length \eqn{s}
containing a sequence of the values 0 or 1. For example, when \eqn{s = 2}, 
the possible states are \eqn{[11]}, \eqn{[10]}, \eqn{[01]}, or \eqn{[00]}, 
corresponding to both species present, only species 1 or species 2 present, 
or both species absent, respectively. The latent state modeled as a 
multivariate Bernoulli random variable:

\deqn{\mathbf{Z}_i \sim \textrm{MVB}(\boldsymbol{\psi}_i)}{Z_i ~ MVB(psi_i)}

where \eqn{\boldsymbol{\psi}_i}{\psi_i} is a vector of length \eqn{2^s} containing the probability of each possible combination of 0s and 1s, such that 
\eqn{\sum\boldsymbol{\psi}_i = 1}{sum(psi_i) = 1}.

For \eqn{s = 2}, the corresponding natural parameters \eqn{f} are

\deqn{f_1 = \log\left(\frac{\psi_{10}}{\psi_{00}}\right)}{f_1 = log(psi_10/psi_00)}
\deqn{f_2 = \log\left(\frac{\psi_{01}}{\psi_{00}}\right)}{f_2 = log(psi_01/psi_00)}
\deqn{f_{12} = \log\left(\frac{\psi_{11}\psi_{00}}{\psi_{10}\psi_{01}}\right)}{f_12 = log((psi_11 * psi_00)/(psi_10 * psi_01))}

The natural parameters can then be modeled as linear functions of covariates.
Covariates for each \eqn{f} must be specified with the \code{stateformulas} argument, 
which takes a character vector of individual formulas of length equal to the number of 
natural parameters (which in turn depends on the number of species in the model).

The observation process is similar to the standard single-species occupancy
model, except that the observations \eqn{\mathbf{y}_{ij}}{y_ij} at site \eqn{i} 
on occasion \eqn{j} are vectors of length \eqn{s} and there are independent
values of detection probability \eqn{p} for each species \eqn{s}:

\deqn{\mathbf{y}_{ij}|\mathbf{Z}_i \sim \textrm{MVB}(\mathbf{Z}_i p_{sij})}{y_ij | Z_i ~ Bernoulli(Z_i  * p_sij)}

Independent detection models (potentially containing different covariates) 
must be provided for each species with the \code{detformulas} argument, 
which takes a character vector of individual formulas with length equal to 
the number of species \eqn{s}.

If you are having problems with separation or boundary estimates (indicated by 
very large parameter estimates and SEs), use of penalized likelihood may
help: see Clipp et al. (2021). \code{occuMulti} supports use of the
Bayes-inspired penalty of Hutchinson et al. (2015). You can set the penalty
value manually using the \code{penalty} argument, or identify the optimal
penalty using K-fold cross validation with the \code{optimizePenalty} function.
See example below.

}

\value{unmarkedFitOccuMulti object describing the model fit.}

\references{

Clipp, H. L., Evans, A., Kessinger, B. E., Kellner, K. F., and C. T. Rota.
  2021. A penalized likelihood for multi-species occupancy models improves 
  predictions of species interactions. Ecology.

Hutchinson, R. A., J. V. Valente, S. C. Emerson, M. G. Betts, and
  T. G. Dietterich. 2015. Penalized Likelihood Methods Improve Parameter
  Estimates in Occupancy Models. Methods in Ecology and Evolution. DOI:
  10.1111/2041-210X.12368
  
MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

Rota, C.T., et al. 2016. A multi-species occupancy model for two or more
  interacting species. Methods in Ecology and Evolution 7: 1164-1173.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuMulti}}}


\examples{

\dontrun{
#Simulate 3 species data
N <- 1000
nspecies <- 3
J <- 5

occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
names(occ_covs) <- paste('occ_cov',1:10,sep='')

det_covs <- list()
for (i in 1:nspecies){
  det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
}
names(det_covs) <- paste('det_cov',1:nspecies,sep='')

#True vals
beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
f1 <- beta[1] + beta[2]*occ_covs$occ_cov1
f2 <- beta[3] + beta[4]*occ_covs$occ_cov2
f3 <- beta[5] + beta[6]*occ_covs$occ_cov3
f4 <- beta[7]
f5 <- beta[8]
f6 <- beta[9]
f7 <- beta[10]
f <- cbind(f1,f2,f3,f4,f5,f6,f7)
z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
colnames(z) <- paste('sp',1:nspecies,sep='')
dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)

psi <- exp(f \%*\% t(dm))
psi <- psi/rowSums(psi)

#True state
ztruth <- matrix(NA,nrow=N,ncol=nspecies)
for (i in 1:N){
  ztruth[i,] <- as.matrix(z[sample(8,1,prob=psi[i,]),])
}

p_true <- c(0.6,0.7,0.5)

# fake y data
y <- list()

for (i in 1:nspecies){
  y[[i]] <- matrix(NA,N,J)
  for (j in 1:N){
    for (k in 1:J){
      y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
    }
  }
}
names(y) <- c('coyote','tiger','bear')

#Create the unmarked data object
data = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)

#Summary of data object
summary(data)
plot(data)

# Look at f parameter design matrix
data@fDesign

# Formulas for state and detection processes

# Length should match number/order of columns in fDesign
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3','~1','~1','~1','~1')

#Length should match number/order of species in data@ylist
detFormulas <- c('~1','~1','~1')

fit <- occuMulti(detFormulas,occFormulas,data)

#Look at output
fit

plot(fit)

#Compare with known values
cbind(c(beta,log(p_true/(1-p_true))),fit@opt$par)

#predict method
lapply(predict(fit,'state'),head)
lapply(predict(fit,'det'),head)

#marginal occupancy
head(predict(fit,'state',species=2))
head(predict(fit,'state',species='bear'))
head(predict(fit,'det',species='coyote'))

#probability of co-occurrence of two or more species
head(predict(fit, 'state', species=c('coyote','tiger')))

#conditional occupancy
head(predict(fit,'state',species=2,cond=3)) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='bear')) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='-bear')) #bear absent
head(predict(fit,'state',species='tiger',cond=c('coyote','-bear')))

#residuals (by species)
lapply(residuals(fit),head)

#ranef (by species)
ranef(fit, species='coyote')

#parametric bootstrap
bt <- parboot(fit,nsim=30)

#update model
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3','~1','~1','~1','~1')
fit2 <- update(fit,stateformulas=occFormulas)

#List of fitted models
fl <- fitList(fit,fit2)
coef(fl)

#Model selection
modSel(fl)

#Fit model while forcing some natural parameters to be 0
#For example: fit model with no species interactions
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3','0','0','0','0')
fit3 <- occuMulti(detFormulas,occFormulas,data)

#Alternatively, you can force all interaction parameters above a certain
#order to be zero with maxOrder. This will be faster.
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3')
fit4 <- occuMulti(detFormulas,occFormulas,data,maxOrder=1)

#Add Bayes penalty term to likelihood. This is useful if your parameter
#estimates are very large, eg because of separation.
fit5 <- occuMulti(detFormulas, occFormulas, data, penalty=1)

#Find optimal penalty term value from a range of possible values using
#K-fold cross validation, and re-fit the model
fit_opt <- optimizePenalty(fit5, penalties=c(0,1,2))
}

}

\keyword{models}
