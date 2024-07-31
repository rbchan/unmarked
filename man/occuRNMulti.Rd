\name{occuRNMulti}

\alias{occuRNMulti}

\title{Fit the Twining et al. (2024) Abundance-mediated Species Interaction Model}

\usage{occuRNMulti(detformulas, stateformulas, data, modelOccupancy,
                   K = rep(10, length(stateformulas)),
                   starts, method="BFGS", se=TRUE, threads = 1, ...)}

\arguments{
    \item{detformulas}{Named list of formulas for the detection models, one per species.}
    \item{stateformulas}{A named list of formulas for the state model. The list should
     have one element per species. Each list element is a single formula (for a dominant species)
     or a list of formulas (for a subordinate species). In the latter case the first formula 
     corresponds to the part of the linear predictor unaffected by interactions,
     and each subsequent formula corresponding to the linear predictor associated
     with interactions with another species. See examples.} 
    \item{data}{An \code{\link{unmarkedFrameOccuMulti}} object.}
    \item{modelOccupancy}{An (optional) character vector of species names,
      corresponding to species for which you want to model occupancy as the
      state parameter instead of abundance. This is only possible for
      subordinate species that do not have an effect on any other species.}
    \item{K}{Vector of upper summation indices (one per species) used to 
      numerically integrate out the latent abundance. This should be set high 
      enough so that it does not affect the parameter estimates. Computation 
      time will increase with K.}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{Logical specifying whether or not to compute standard
      errors.}
    \item{threads}{Set the number of threads to use for optimization in C++, if
      OpenMP is available on your system. Increasing the number of threads
      may speed up optimization in some cases by running the likelihood 
      calculation in parallel. If \code{threads=1} (the default), OpenMP is disabled.} 
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{Fit the Twining et al. (2024) abundance-mediated species interaction model
  using multi-species detection/non-detection data. The model uses the Royle-Nichols
  model (Royle and Nichols 2003) to estimate abundance from detection/non-detection data.
  The abundance of a dominant species can affect the abundance of subordinate
  species, and this interaction can be mediated by covariates. A maximum
  of 3 species are allowed.}
 
\details{

See \code{\link{unmarkedFrame}} and \code{\link{unmarkedFrameOccuMulti}} for a
description of how to supply data to the \code{data} argument.

Royle and Nichols (2003) described how species abundance at a site can be 
obtained from detection/non-detection data by exploiting the expected positive correlation
between number of individuals of a species present at a site \eqn{i} (\eqn{N_i}) 
and overall detection probability of at least one individual (\eqn{p_i}):

\deqn{p_i = 1 - (1 - r_i)^{N_i}}

where \eqn{r_i} is the probability of detecting an individual animal. Abundance
is then typically assumed to come from a Poisson distribution with parameter
\eqn{lambda_i} which can be a function of covariates.

Twining et al. (2024) extend this model to multi-species data and additionally
model interaction(s) between a dominant species (\eqn{D}) and a subordinate 
species (\eqn{S}) by including parameter(s) in the linear predictor for \eqn{\lambda_S}
that are interacted with \eqn{N_D}. For example, expected abundance of the subordinate
species can include an interaction that is dependent on some covariate \eqn{x}:

\deqn{\mathrm{log}(\lambda_{S,i}) = \beta_0 + (\gamma_0 + \gamma_x \cdot x_i) \cdot N_{D,i}}

where the \eqn{gamma} terms are the interaction parameters. Note that the 
subordinate species \eqn{S} could not itself affect abundance of dominant species
\eqn{D}; thus the interactions are directional and not symmetric as with the
Rota et al. (2016) model. Subordinate species could have occupancy models instead of 
abundance models (controlled with function argument \code{modelOccupancy}),
in which case the parameter of interest is \eqn{\psi_{S}}. 
However, any species for which you want to model an effect on another species 
(i.e., a dominant species), you cannot fit an occupancy model.

In principle this model can work with any number of of species and interactions.
\code{occuRNMulti} handles only a maximum of 3 species. However,
it supports all possible interaction combinations for 2 or 3 species. 
The interactions are specified using the \code{stateformulas} argument.
See the examples for details.

Detection models for each species are independent and potentially include different
covariates. Detection models are provided to the \code{detformulas} argument, 
which takes a named list of individual formulas with length equal to 
the number of species.
}

\value{unmarkedFitOccuRNMulti object describing the model fit.}

\references{
  Twining 2024
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarkedFrameOccuMulti}}, \code{\link{occuRN}}}


\examples{

\donttest{

# A model where a dominant species (coyote) affects an intermediate species 
# (marten) which affects a subordinate species (fisher)
# e.g., sp1 --> sp2 --> sp3
set.seed(1)

M <- 200
J <- 5

N1 <- N2 <- N3 <- 25
while((max(N1) > 20) | (max(N2) > 20) | (max(N3) > 20)){
# Covariates
sc <- data.frame(x1 = rnorm(M), x2 = rnorm(M))
oc <- data.frame(x3 = rnorm(M*J))

# Coyote expected abundance
l1 <- 2

# True coyote abundance
N1 <- rpois(M, l1)

# Lambda for marten
l2 <- exp(log(1.3) - 0.3*N1 + 0.5*sc$x1*N1)

# True marten abundance
N2 <- rpois(M, l2)

# Lambda for fisher
l3 <- exp(log(1) + -0.2 * N2)
N3 <- rpois(M, l3)
}

# True detection prob for individuals of both species
rs <- 0.2
ps1 <- 1 - (1-rs)^N1
ps2 <- 1 - (1-rs)^N2
ps3 <- 1 - (1-rs)^N3

# True parameter values
truth <- c(log(2), log(1.3), -0.3, 0.5, log(1), -0.2, rep(log(0.2/(1-0.2)), 3)) 

# Simulate observations
y1 <- y2 <- y3 <- matrix(NA, M, J)
for (m in 1:M){
  for (j in 1:J){
    y1[m,] <- rbinom(J, 1, ps1[m])
    y2[m,] <- rbinom(J, 1, ps2[m])
    y3[m,] <- rbinom(J, 1, ps3[m])
  }
}

umf <- unmarkedFrameOccuMulti(y=list(coyote=y1, marten=y2, fisher=y3),
                              siteCovs=sc, obsCovs=oc)

# Specify formulas for model as named list
# sp1 --> sp2 --> sp3
# Coyote is unaffected by other species, so it includes just a single formula
# Marten is affected by coyote, so it has a list where the first formula
# is the part of the linear predictor not affected by coyote (here intercept-only),
# and the second (named) list element is the coyote effect, which is mediated
# by covariate x1.
# Finally, fisher has an intercept-only model for the independent part,
# and is affected by marten abundance, but without any covariates.
sf <- list(coyote = ~1,
           marten = list(~1, coyote = ~x1),
           fisher = list(~1, marten = ~1))

df <- list(coyote = ~1, marten = ~1, fisher=~1)

K <- c(max(N1) + 3, max(N2) + 3, max(N3) + 3)

fit <- occuRNMulti(df, sf, umf, K=K, threads=4)
summary(fit)

cbind(unm=coef(fit), truth=truth)

# Predicted abundance for each species at each site
# This plugs in expected abundance (lambda) of dominant species
# into the linear predictor, rather than using the latent abundance N
pr <- predict(fit, type='state')
lapply(pr, head)

# Detection probability for each species
pr <- predict(fit, type = "state", level=NULL)

# Example formula structures for other models
# (could add covariates to any formula)
# sp1 --> sp2
sf <- list(coyote = ~1, marten = list(~1, coyote = ~1))
# sp2 <-- sp1 --> sp3
sf <- list(coyote = ~1,
           marten = list(~1, coyote = ~1),
           fisher = list(~1, coyote = ~1))
# sp1 --> sp3 <-- sp2
sf <- list(coyote = ~1,
           marten = ~1,
           fisher = list(~1, coyote = ~1, marten = ~1))
# sp1 --> sp2 --> sp3 AND sp1 --> sp3
sf <- list(coyote = ~1,
           marten = list(~1, coyote = ~1),
           fisher = list(~1, coyote = ~1, marten = ~1))
}

}
