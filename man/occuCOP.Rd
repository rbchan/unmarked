\name{occuCOP}

\alias{occuCOP}

\encoding{UTF-8}

\title{Fit the occupancy model using count dta}

\usage{
occuCOP(data, 
        psiformula = ~1, lambdaformula = ~1, 
        psistarts, lambdastarts, starts,
        method = "BFGS", se = TRUE, 
        engine = c("C", "R"), na.rm = TRUE, 
        return.negloglik = NULL, L1 = FALSE, ...)}

\arguments{

    \item{data}{An \code{\link{unmarkedFrameOccuCOP}} object created with the \code{\link{unmarkedFrameOccuCOP}} function.}
    
    \item{psiformula}{Formula describing the occupancy covariates.}
    
    \item{lambdaformula}{Formula describing the detection covariates.}
    
    \item{psistarts}{Vector of starting values for likelihood maximisation with \code{\link{optim}} for occupancy probability \eqn{\psi}{psi}. These values must be logit-transformed (with \code{\link{qlogis}}) (see details). By default, optimisation will start at 0, corresponding to an occupancy probability of 0.5 (\code{plogis(0)} is 0.5).}
    
    \item{lambdastarts}{Vector of starting values for likelihood maximisation with \code{\link{optim}} for detection rate \eqn{\lambda}{lambda}. These values must be log-transformed (with \code{\link{log}}) (see details). By default, optimisation will start at 0, corresponding to detection rate of 1 (\code{exp(0)} is 1).}
    
    \item{starts}{Vector of starting values for likelihood maximisation with \code{\link{optim}}. If \code{psistarts} and \code{lambdastarts} are provided, \code{starts = c(psistarts, lambdastarts)}.}

    \item{method}{Optimisation method used by \code{\link{optim}}.}
    
    \item{se}{Logical specifying whether to compute (\code{se=TRUE}) standard errors or not (\code{se=FALSE}).}
    
    \item{engine}{Code to use for optimisation. Either \code{"C"} for fast C++ code, or \code{"R"} for native R code.}
    
    \item{na.rm}{Logical specifying whether to fit the model (\code{na.rm=TRUE}) or not (\code{na.rm=FALSE}) if there are NAs in the \code{\link{unmarkedFrameOccuCOP}} object.}
    
    \item{return.negloglik}{A list of vectors of parameters (\code{c(psiparams, lambdaparams)}). If specified, the function will not maximise likelihood but return the negative log-likelihood for the those parameters in the \code{nll} column of a dataframe. See an example below.}

    \item{L1}{Logical specifying whether the length of observations (\code{L}) are purposefully set to 1 (\code{L1=TRUE}) or not (\code{L1=FALSE}).}

    \item{\dots}{Additional arguments to pass to \code{\link{optim}}, such as lower and upper bounds or a list of control parameters.}
  }

\description{This function fits a single season occupancy model using count data.}

\details{
  
  See \code{\link{unmarkedFrameOccuCOP}} for a description of how to supply data to the \code{data} argument. See \code{\link{unmarkedFrame}} for a more general documentation of \code{unmarkedFrame} objects for the different models implemented in \pkg{unmarked}.
  
  \subsection{The COP occupancy model}{
    
    \code{occuCOP} fits a single season occupancy model using count data, as described in Pautrel et al. (2023). 

    The \strong{occupancy sub-model} is:
    
    \deqn{z_i \sim \text{Bernoulli}(\psi_i)}{z_i ~ Bernoulli(psi_i)}
    
    \itemize{
      \item With \eqn{z_i}{z_i} the occupany state of site \eqn{i}{i}. \eqn{z_i=1}{z_i = 1} if site \eqn{i}{i} is occupied by the species, \emph{i.e.} if the species is present in site \eqn{i}{i}. \eqn{z_i=0}{z_i = 0} if site \eqn{i}{i} is not occupied.
      \item With \eqn{\psi_i}{psi_i} the occupancy probability of site \eqn{i}{i}.
    }
    
    The \strong{observation sub-model} is:
    
    \deqn{
      N_{ij} | z_i = 1 \sim \text{Poisson}(\lambda_{ij} L_{ij}) \\
      N_{ij} | z_i = 0 \sim 0
    }{
      N_ij | z_i = 1 ~ Poisson(lambda_is*L_is)
      N_ij | z_i = 0 ~ 0
    }
    
    \itemize{
      \item With \eqn{N_{ij}}{N_ij} the count of detection events in site \eqn{i}{i} during observation \eqn{j}{j}.
      \item With \eqn{\lambda_{ij}}{lambda_ij} the detection rate in site \eqn{i}{i} during observation \eqn{j}{j} (\emph{for example, 1 detection per day.}).
      \item With \eqn{L_{ij}}{L_ij} the length of observation \eqn{j}{j} in site \eqn{i}{i} (\emph{for example, 7 days.}).
    }
    
    What we call "observation" (\eqn{j}{j}) here can be a sampling occasion, a transect, a discretised session. Consequently, the unit of \eqn{\lambda_{ij}}{lambda_ij} and \eqn{L_{ij}}{L_ij} can be either a time-unit (day, hour, ...) or a space-unit (kilometer, meter, ...).
  }
  
  \subsection{The transformation of parameters \eqn{\psi} and \eqn{\lambda}}{
    In order to perform unconstrained optimisation, parameters are transformed.
    
    The occupancy probability (\eqn{\psi}) is transformed with the logit function (\code{psi_transformed = qlogis(psi)}). It can be back-transformed with the "inverse logit" function (\code{psi = plogis(psi_transformed)}).

    The detection rate (\eqn{\lambda}) is transformed with the log function (\code{lambda_transformed = log(lambda)}). It can be back-transformed with the exponential function (\code{lambda = exp(lambda_transformed)}).
  }
  
}

\value{\code{unmarkedFitOccuCOP} object describing the model fit. See the \code{\linkS4class{unmarkedFit}} classes.}

\references{

Pautrel, L., Moulherat, S., Gimenez, O. & Etienne, M.-P. Submitted. \emph{Analysing biodiversity observation data collected in continuous time: Should we use discrete or continuous-time occupancy models?} Preprint at \doi{10.1101/2023.11.17.567350}.

}

\author{Léa Pautrel}

\seealso{
  \code{\link{unmarked}}, 
  \code{\link{unmarkedFrameOccuCOP}},
  \code{\link{unmarkedFit-class}}
}


\examples{
set.seed(123)
options(max.print = 50)

# We simulate data in 100 sites with 3 observations of 7 days per site.
nSites <- 100
nObs <- 3

# For an occupancy covariate, we associate each site to a land-use category.
landuse <- sample(factor(c("Forest", "Grassland", "City"), ordered = TRUE), 
                  size = nSites, replace = TRUE)
simul_psi <- ifelse(landuse == "Forest", 0.8, 
                    ifelse(landuse == "Grassland", 0.4, 0.1))
z <- rbinom(n = nSites, size = 1, prob = simul_psi)

# For a detection covariate, we create a fake wind variable.
wind <- matrix(rexp(n = nSites * nObs), nrow = nSites, ncol = nObs)
simul_lambda <- wind / 5
L = matrix(7, nrow = nSites, ncol = nObs)

# We now simulate count detection data
y <- matrix(rpois(n = nSites * nObs, lambda = simul_lambda * L), 
            nrow = nSites, ncol = nObs) * z

# We create our unmarkedFrameOccuCOP object
umf <- unmarkedFrameOccuCOP(
  y = y,
  L = L,
  siteCovs = data.frame("landuse" = landuse),
  obsCovs = list("wind" = wind)
)
print(umf)

# We fit our model without covariates
fitNull <- occuCOP(data = umf)
print(fitNull)

# We fit our model with covariates
fitCov <- occuCOP(data = umf, psiformula = ~ landuse, lambdaformula = ~ wind)
print(fitCov)

# We back-transform the parameter's estimates
## Back-transformed occupancy probability with no covariates
backTransform(fitNull, "psi")

## Back-transformed occupancy probability depending on habitat use
predict(fitCov,
        "psi",
        newdata = data.frame("landuse" = c("Forest", "Grassland", "City")),
        appendData = TRUE)

## Back-transformed detection rate with no covariates
backTransform(fitNull, "lambda")

## Back-transformed detection rate depending on wind
predict(fitCov,
        "lambda",
        appendData = TRUE)

## This is not easily readable. We can show the results in a clearer way, by:
##  - adding the site and observation
##  - printing only the wind covariate used to get the predicted lambda
cbind(
  data.frame(
    "site" = rep(1:nSites, each = nObs),
    "observation" = rep(1:nObs, times = nSites),
    "wind" = getData(fitCov)@obsCovs
  ),
  predict(fitCov, "lambda", appendData = FALSE)
)

# We can choose the initial parameters when fitting our model.
# For psi, intituively, the initial value can be the proportion of sites 
#          in which we have observations.
(psi_init <- mean(rowSums(y) > 0))

# For lambda, the initial value can be the mean count of detection events 
#             in sites in which there was at least one observation.
(lambda_init <- mean(y[rowSums(y) > 0, ]))

# We have to transform them.
occuCOP(
  data = umf,
  psiformula = ~ 1,
  lambdaformula = ~ 1,
  psistarts = qlogis(psi_init),
  lambdastarts = log(lambda_init)
)

# If we have covariates, we need to have the right length for the start vectors.
# psi ~ landuse --> 3 param to estimate: Intercept, landuseForest, landuseGrassland
# lambda ~ wind --> 2 param to estimate: Intercept, wind
occuCOP(
  data = umf,
  psiformula = ~ landuse,
  lambdaformula = ~ wind,
  psistarts = rep(qlogis(psi_init), 3),
  lambdastarts = rep(log(lambda_init), 2)
)

# And with covariates, we could have chosen better initial values, such as the
# proportion of sites in which we have observations per land-use category.
(psi_init_covs <- c(
  "City" = mean(rowSums(y[landuse == "City", ]) > 0),
  "Forest" = mean(rowSums(y[landuse == "Forest", ]) > 0),
  "Grassland" = mean(rowSums(y[landuse == "Grassland", ]) > 0)
))
occuCOP(
  data = umf,
  psiformula = ~ landuse,
  lambdaformula = ~ wind,
  psistarts = qlogis(psi_init_covs))

# We can fit our model with a different optimisation algorithm.
occuCOP(data = umf, method = "Nelder-Mead")

# We can run our model with a C++ or with a R likelihood function.
## They give the same result. 
occuCOP(data = umf, engine = "C", psistarts = 0, lambdastarts = 0)
occuCOP(data = umf, engine = "R", psistarts = 0, lambdastarts = 0)

## The C++ (the default) is faster.
system.time(occuCOP(data = umf, engine = "C", psistarts = 0, lambdastarts = 0))
system.time(occuCOP(data = umf, engine = "R", psistarts = 0, lambdastarts = 0))

## However, if you want to understand how the likelihood is calculated,
## you can easily access the R likelihood function.
print(occuCOP(data = umf, engine = "R", psistarts = 0, lambdastarts = 0)@nllFun)

# Finally, if you do not want to fit your model but only get the likelihood,
# you can get the negative log-likelihood for a given set of parameters.
occuCOP(data = umf, return.negloglik = list(
  c("psi" = qlogis(0.25), "lambda" = log(2)),
  c("psi" = qlogis(0.5), "lambda" = log(1)),
  c("psi" = qlogis(0.75), "lambda" = log(0.5))
))
}

\keyword{models}
