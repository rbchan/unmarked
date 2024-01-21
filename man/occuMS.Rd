\name{occuMS}

\alias{occuMS}

\title{Fit Single-Season and Dynamic Multi-State Occupancy Models}

\usage{occuMS(detformulas, psiformulas, phiformulas=NULL, data, 
    parameterization=c("multinomial","condbinom"),
    starts, method="BFGS", se=TRUE, engine=c("C","R"), silent=FALSE, ...)}

\arguments{
    \item{detformulas}{Character vector of formulas for detection probabilities.
      See details for a description of how to order these formulas.}
    \item{psiformulas}{Character vector of formulas for occupancy probabilities. 
      See details for a description of how to order these formulas.}
    \item{phiformulas}{Character vector of formulas for state transition probabilities. 
      Only used if you are fitting a dynamic model. See details for a 
      description of how to order these formulas.}
    \item{data}{An \code{\link{unmarkedFrameOccuMS}} object}
    \item{parameterization}{Either \code{"multinomial"} for the multinomial 
      parameterization (MacKenzie et al. 2009) which allows an arbitrary
      number of occupancy states, or \code{"condbinom"} for the conditional 
      binomial parameterization (Nichols et al. 2007) which requires exactly
      3 occupancy states. See details.}
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

\description{This function fits single-season and dynamic multi-state occupancy models with both the multinomial and conditional binomial parameterizations.}

\details{

Traditional occupancy models fit data with exactly two states: detection and
non-detection (MacKenzie et al. 2002).
The \code{occuMS} function fits models to occupancy data for which there are 
greater than 2 states (Nichols et al 2007, MacKenzie et al. 2009). For example, 
detections may be further divided into multiple biologically relevant categories, 
e.g. breeding vs. non-breeding, or some/many individuals present. As with
detection status, classification of these additional occupancy states is likely
to be imperfect.

Multiple parameterizations for multi-state occupancy models have been proposed.
The \code{occuMS} function fits two at present: the "conditional binomial" 
parameterization of Nichols et al. (2007), and the more general "multinomial"
parameterization of MacKenzie et al. (2009). Both single-season
and dynamic models are possible with \code{occuMS} (MacKenzie et al. 2009).

The conditional binomial parameterization 
(\code{parameterization = 'condbinom'}) models occupancy and the presence or
absence of an additional biological state of interest given the species
is present (typically breeding status). Thus, there should be exactly 3 occupancy 
states in the data: 0 (non-detection); 1 (detection, no evidence of breeding);
or 2 (detection, evidence of breeding). 

Two state parameters are estimated:
\eqn{\psi}, the probability of occupancy, and \eqn{R}, the probability of 
successful reproduction given an occupied state (although this could be some
other binary biological condition). Covariates (in \code{siteCovs}) can be 
supplied for either or both of these parameters with the \code{stateformulas} 
argument, which takes a character vector of R-style formulas with length = 2, 
with formulas in the order (\eqn{\psi}, \eqn{R}). For example, to fit a model 
where \eqn{\psi} varies with a landcover covariate and \eqn{R} is constant,
\code{stateformulas = c('~landcover','~1')}. 

There are three detection parameters associated with the
conditional binomial parameterization: \eqn{p_1}, the probability of 
detecting the species given true state 1; \eqn{p_2}, the probability of detecting
the species given true state 2; and \eqn{\delta}, the probability of detecting 
state 2 (i.e., breeding), given that the species has been detected.
See MacKenzie et al. (2009), pages 825-826 for more details.
As with occupancy, covariates (in \code{obsCovs}) can be supplied for these 
detection probabilities with the \code{detformulas} argument, which takes a 
character vector of formulas with length = 3 in the order 
(\eqn{p_1}, \eqn{p_2}, \eqn{\delta}). So, to fit a model where \eqn{p_1} varies 
with temperature and the other two parameters are constant, 
\code{detformulas = c('~temp','~1','~1')}.

The multinomial parameterization (\code{parameterization = "multinomial"}) is
more general, allowing an arbitrary number of occupancy states \eqn{S}.
\eqn{S} - 1 occupancy probabilities \eqn{\psi} are estimated. Thus, if there
are \eqn{S} = 4 occupancy states (0, 1, 2, 3), \code{occuMS} estimates \eqn{\psi_1},
\eqn{\psi_2}, and \eqn{\psi_3} (the probability of state 0 can be obtained by 
subtracting the others from 1). Covariates can be supplied for each occupancy
probability with a character vector with length \eqn{S-1}, e.g. 
\code{stateformulas =} \code{c('~landcover','~1','~1')} where \eqn{\psi_1} varies with
landcover and \eqn{\psi_2} and \eqn{\psi_3} are constant.

The number of detection probabilities estimated quickly expands as \eqn{S}
increases, equal to \eqn{S \times (S-1) / 2}. In the simplest case
(when \eqn{S} = 3), there are 3 detection probabilities: \eqn{p_{11}}, 
the probability of detecting state 1 given true state 1; \eqn{p_{12}}, 
the probability of detecting state 1 given true state 2; and \eqn{p_{22}}, 
the probability of detecting state 2 given true state 2. 
Covariates can be supplied for any or all of these detection probabilities with 
the \code{detformulas} argument, which takes a character vector of formulas 
with length = 3 in the order (\eqn{p_{11}}, \eqn{p_{12}}, \eqn{p_{22}}). So, 
to fit a model where \eqn{p_{11}} varies with temperature and the other two detection 
probabilities are constant, \code{detformulas = c('~temp','~1','~1')}.
If there were \eqn{S} = 4 occupancy states, there are 6 estimated detection
probabilities and the order is (\eqn{p_{11}}, \eqn{p_{12}}, \eqn{p_{13}},
\eqn{p_{22}}, \eqn{p_{23}}, \eqn{p_{33}}), and so on. See MacKenzie et al. (2009)
for a more detailed explanation.

Dynamic (multi-season) models can be fit as well for both parameterizations
(MacKenzie et al. 2009). In a standard dynamic occupancy model, additional
parameters for probabilities of colonization (i.e., state 0 -> 1) and 
extinction (1 -> 0) are estimated. In a multi-state context, we must estimate a
transition probability matrix (\eqn{\phi}) between all possible states. You can 
provide formulas for some of the probabilities in this matrix using the 
\code{phiformulas} argument. The approach differs depending on parameterization.

For the conditional binomial parameterization, \code{phiformulas} is a 
character vector of length 6. The first three elements are formulas for the 
probability a site is occupied at time \eqn{t} given that it was previously
in states 0, 1, or 2 at time \eqn{t-1} (\code{phi0, phi1, phi2}). Elements 4-6 
are formulas for the probability of reproduction (or other biological state) 
given state 0, 1, or 2 at time \eqn{t-1} (\code{R0, R1, R2}). See 
\code{umf@phiOrder$cond_binom} for a reminder of the correct order, where
\code{umf} is your \code{unmarkedFrameOccuMS}.

For the multinomial parameterization, \code{phiformulas} can be used to provide
formulas for some transitions between different occupancy states. You can't 
give formulas for the probabilities of remaining in the same state between 
seasons to keep the model identifiable. Thus, if there are 3 possible states
(0, 1, 2), \code{phiformulas} should contain 6 formulas for the following 
transitions: \code{p(0->1), p(0->2), p(1->0), p(1->2), p(2->0), p(2->1)},
in that order (and similar for more than 3 states). The remaining probabilities
of staying in the same state between seasons can be obtained via subtraction.
See \code{umf@phiOrder$multinomial} for the correct order matching the number
of states in your dataset.

See \code{\link{unmarkedFrameOccuMS}} for a description of how to supply data 
to the \code{data} argument.
}

\value{unmarkedFitOccuMS object describing the model fit.}

\references{
MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

MacKenzie, D. I., Nichols, J. D., Seamans, M. E., and R. J. Gutierrez, 2009. 
  Modeling species occurrence dynamics with multiple states and imperfect 
  detection. Ecology 90: 823-835.

Nichols, J. D., Hines, J. E., Mackenzie, D. I., Seamans, M. E., and 
  R. J. Gutierrez. 2007. Occupancy estimation and modeling with multiple states 
  and state uncertainty. Ecology 88: 1395-1400.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuMS}}}


\examples{

\dontrun{

#Simulate data

#Parameters
N <- 500; J <- 5; S <- 3
site_covs <- matrix(rnorm(N*2),ncol=2)
obs_covs <- matrix(rnorm(N*J*2),ncol=2)
a1 <- -0.5; b1 <- 1; a2 <- -0.6; b2 <- -0.7

##################################
## Multinomial parameterization ##
##################################

p11 <- -0.4; p12 <- -1.09; p22 <- -0.84
truth <- c(a1,b1,a2,b2,p11,0,p12,p22)

#State process
lp <- matrix(NA,ncol=S,nrow=N)
for (n in 1:N){
  lp[n,2] <- exp(a1+b1*site_covs[n,1])
  lp[n,3] <- exp(a2+b2*site_covs[n,2])
  lp[n,1] <- 1  
}
psi_mat <- lp/rowSums(lp)

z <- rep(NA,N)
for (n in 1:N){
  z[n] <- sample(0:2, 1, replace=T, prob=psi_mat[n,])
}

probs_raw <- matrix(c(1,0,0,1,exp(p11),0,1,exp(p12),exp(p22)),nrow=3,byrow=T)
probs_raw <- probs_raw/rowSums(probs_raw)
  
y <- matrix(0,nrow=N,ncol=J)
for (n in 1:N){

  probs <- switch(z[n]+1,
                  probs_raw[1,],
                  probs_raw[2,],
                  probs_raw[3,])
  if(z[n]>0){
    y[n,] <- sample(0:2, J, replace=T, probs)
  }
}

#Construct unmarkedFrame
umf <- unmarkedFrameOccuMS(y=y,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs))

#Formulas

#3 states, so detformulas is a character vector of formulas of 
#length 3 in following order:
#1) p[11]: prob of detecting state 1 given true state 1
#2) p[12]: prob of detecting state 1 given true state 2
#3) p[22]: prob of detecting state 2 given true state 2
detformulas <- c('~V1','~1','~1')
#If you had 4 states, it would be p[11],p[12],p[13],p[22],p[23],p[33] and so on

#3 states, so stateformulas is a character vector of length 2 in following order:
#1) psi[1]: probability of state 1
#2) psi[2]: probability of state 2
#You can get probability of state 0 (unoccupied) as 1 - psi[1] - psi[2]
stateformulas <- c('~V1','~V2')

#Fit model
fit <- occuMS(detformulas, stateformulas, data=umf,
              parameterization="multinomial")

#Look at results
fit
#Compare with truth
cbind(truth=truth,estimate=coef(fit))

#Generate predicted values
lapply(predict(fit,type='psi'),head)
lapply(predict(fit,type='det'),head)

#Fit a null model
detformulas <- rep('~1',3)
stateformulas <- rep('~1',2)
fit_null <- occuMS(detformulas, stateformulas, data=umf,
                   parameterization="multinomial")

#Compare fits
modSel(fitList(fit,fit_null))

###########################################
## Conditional binomial parameterization ##
###########################################

p11 <- 0.4; p12 <- 0.6; p22 <- 0.8
truth_cb <- c(a1,b1,a2,b2,qlogis(p11),0,qlogis(c(p12,p22)))

#Simulate data

#State process
psi_mat <- matrix(NA,ncol=S,nrow=N)
for (n in 1:N){
  psi_mat[n,2] <- plogis(a1+b1*site_covs[n,1])
  psi_mat[n,3] <- plogis(a2+b2*site_covs[n,2])
}
psi_bin <- matrix(NA,nrow=nrow(psi_mat),ncol=ncol(psi_mat))
psi_bin[,1] <- 1-psi_mat[,2]
psi_bin[,2] <- (1-psi_mat[,3])*psi_mat[,2]
psi_bin[,3] <- psi_mat[,2]*psi_mat[,3]
z <- rep(NA,N)
for (n in 1:N){
  z[n] <- sample(0:2, 1, replace=T, prob=psi_bin[n,])
}

#Detection process
y_cb <- matrix(0,nrow=N,ncol=J)
for (n in 1:N){
  #p11 = p1; p12 = p2; p22 = delta
  probs <- switch(z[n]+1,
                  c(1,0,0),
                  c(1-p11,p11,0),
                  c(1-p12,p12*(1-p22),p12*p22)) 
  if(z[n]>0){
    y_cb[n,] <- sample(0:2, J, replace=T, probs)
  }
}

#Build unmarked frame
umf2 <- unmarkedFrameOccuMS(y=y_cb,siteCovs=as.data.frame(site_covs),
                           obsCovs=as.data.frame(obs_covs))

#Formulas

#detformulas is a character vector of formulas of length 3 in following order:
#1) p[1]: prob of detecting species given true state 1
#2) p[2]: prob of detecting species given true state 2
#3) delta: prob of detecting state 2 (eg breeding) given species was detected
detformulas <- c('~V1','~1','~1')

#stateformulas is a character vector of length 2 in following order:
#1) psi: probability of occupancy
#2) R: probability state 2 (eg breeding) given occupancyc
stateformulas <- c('~V1','~V2')

#Fit model
fit_cb <- occuMS(detformulas, stateformulas, data=umf2,
                 parameterization='condbinom')

#Look at results
fit_cb
#Compare with truth
cbind(truth=truth_cb,estimate=coef(fit_cb))

#Generate predicted values
lapply(predict(fit_cb,type='psi'),head)
lapply(predict(fit_cb,type='det'),head)


##################################
## Dynamic (multi-season) model ##
##################################

#Simulate data-----------------------------------------------
N <- 500 #Number of sites
T <- 3 #Number of primary periods
J <- 5 #Number of secondary periods
S <- 3 #Number of occupancy states (0,1,2)

#Generate covariates
site_covs <- as.data.frame(matrix(rnorm(N*2),ncol=2))
yearly_site_covs <- as.data.frame(matrix(rnorm(N*T*2),ncol=2))
obs_covs <- as.data.frame(matrix(rnorm(N*J*T*2),ncol=2))

#True parameter values
b <- c(
  #Occupancy parameters
  a1=-0.5, b1=1, a2=-0.6, b2=-0.7,
  #Transition prob (phi) parameters
  phi01=0.7, phi01_cov=-0.5, phi02=-0.5, phi10=1.2, 
  phi12=0.3, phi12_cov=1.1, phi20=-0.3, phi21=1.4, phi21_cov=0,
  #Detection prob parameters
  p11=-0.4, p11_cov=0, p12=-1.09, p22=-0.84
)

#Generate occupancy probs (multinomial parameterization)
lp <- matrix(1, ncol=S, nrow=N)
lp[,2] <- exp(b[1]+b[2]*site_covs[,1])
lp[,3] <- exp(b[3]+b[4]*site_covs[,2])
psi <- lp/rowSums(lp)

#True occupancy state matrix
z <- matrix(NA, nrow=N, ncol=T)

#Initial occupancy
for (n in 1:N){
  z[n,1] <- sample(0:(S-1), 1, prob=psi[n,])
}

#Raw phi probs
phi_raw <- matrix(NA, nrow=N*T, ncol=S^2-S)
phi_raw[,1] <- exp(b[5]+b[6]*yearly_site_covs[,1]) #p[0->1]
phi_raw[,2] <- exp(b[7]) #p[0->2]
phi_raw[,3] <- exp(b[8]) #p[1->0]
phi_raw[,4] <- exp(b[9]+b[10]*yearly_site_covs[,2]) #p[1->2]
phi_raw[,5] <- exp(b[11]) #p[2->0]
phi_raw[,6] <- exp(b[12]+b[13]*yearly_site_covs[,1])

#Generate states in times 2..T
px <- 1
for (n in 1:N){
  for (t in 2:T){
    phi_mat <- matrix(c(1, phi_raw[px,1], phi_raw[px,2],  # phi|z=0
                        phi_raw[px,3], 1, phi_raw[px,4],  # phi|z=1
                        phi_raw[px,5], phi_raw[px,6], 1), # phi|z=2
                      nrow=S, byrow=T)
    phi_mat <- phi_mat/rowSums(phi_mat)
    z[n, t] <- sample(0:(S-1), 1, prob=phi_mat[z[n,(t-1)]+1,])
    px <- px + 1
    if(t==T) px <- px + 1 #skip last datapoint for each site
  }
}

#Raw p probs
p_mat <- matrix(c(1, 0, 0, #p|z=0
                  1, exp(b[14]), 0, #p|z=1
                  1, exp(b[16]), exp(b[17])), #p|z=2 
                nrow=S, byrow=T)
p_mat <- p_mat/rowSums(p_mat)

#Simulate observation data
y <- matrix(0, nrow=N, ncol=J*T)
for (n in 1:N){
  yx <- 1
  for (t in 1:T){
    if(z[n,t]==0){
      yx <- yx + J
      next
    }
    for (j in 1:J){
      y[n, yx] <- sample(0:(S-1), 1, prob=p_mat[z[n,t]+1,])
      yx <- yx+1
    }
  }
}
#-----------------------------------------------------------------

#Model fitting

#Build UMF
umf <- unmarkedFrameOccuMS(y=y, siteCovs=site_covs,
                           obsCovs=obs_covs,
                           yearlySiteCovs=yearly_site_covs,
                           numPrimary=3)
summary(umf)

#Formulas
#Initial occupancy
psiformulas <- c('~V1','~V2') #on psi[1] and psi[2]

#Transition probs
#Guide to order:
umf@phiOrder$multinomial
phiformulas <- c('~V1','~1','~1','~V2','~1','~V1')

#Detection probability
detformulas <- c('~V1','~1','~1') #on p[1|1], p[1|2], p[2|2]

#Fit model
(fit <- occuMS(detformulas=detformulas, psiformulas=psiformulas,
              phiformulas=phiformulas, data=umf))

#Compare with truth
compare <- cbind(b,coef(fit),
                 coef(fit)-1.96*SE(fit),coef(fit)+1.96*SE(fit))
colnames(compare) <- c('truth','estimate','lower','upper')
round(compare,3)

#Estimated phi matrix for site 1
phi_est <- predict(fit, 'phi', se.fit=F)
phi_est <- sapply(phi_est, function(x) x$Predicted[1])
phi_est_mat <- matrix(NA, nrow=S, ncol=S)
phi_est_mat[c(4,7,2,8,3,6)] <- phi_est
diag(phi_est_mat) <- 1 - rowSums(phi_est_mat,na.rm=T)

#Actual phi matrix for site 1
phi_act_mat <- diag(S)
phi_act_mat[c(4,7,2,8,3,6)] <- phi_raw[1,]
phi_act_mat <- phi_act_mat/rowSums(phi_act_mat)

#Compare
cat('Estimated phi\n')
phi_est_mat
cat('Actual phi\n')
phi_act_mat

#Rough check of model fit
fit_sim <- simulate(fit, nsim=20)
hist(sapply(fit_sim,mean),col='gray')
abline(v=mean(umf@y),col='red',lwd=2)
#line should fall near middle of histogram

}

}

\keyword{models}
