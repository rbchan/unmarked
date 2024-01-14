\name{predict}
\docType{methods}
\alias{predict}
\alias{predict-methods}
\alias{predict,unmarkedFit-method}
\alias{predict,unmarkedFitOccuMulti-method}
\alias{predict,unmarkedFitOccuMS-method}
\alias{predict,unmarkedFitList-method}
\alias{predict,unmarkedRanef-method}

\title{Methods for Function predict in Package 'unmarked'}

\description{
These methods return predicted values from \code{unmarkedFit} objects, \code{fitList}s,
or \code{ranef} output. Most object types use the same method, but some have 
unique options; see usage below.
}

\usage{
\S4method{predict}{unmarkedFit}(object, type, newdata, backTransform = TRUE,
  na.rm = TRUE, appendData = FALSE, level = 0.95, re.form = NULL, ...)
\S4method{predict}{unmarkedFitOccuMulti}(object, type, newdata, se.fit = TRUE,
  level = 0.95, species = NULL, cond = NULL, nsims = 100, ...)
\S4method{predict}{unmarkedFitOccuMS}(object, type, newdata, se.fit = TRUE,
  level = 0.95, nsims = 100, ...)
\S4method{predict}{unmarkedFitList}(object, type, newdata = NULL,
  backTransform = TRUE, appendData = FALSE, level = 0.95)
\S4method{predict}{unmarkedRanef}(object, func, nsims = 100, ...)
}

\arguments{
  \item{object}{A \code{unmarkedFit}, \code{unmarkedFitList}, or \code{unmarkedRanef} object.}
  \item{type}{The submodel to predict values for, such as \code{state} or \code{det}.
    The available types depend on the model.}
  \item{newdata}{Optional; provide a \code{data.frame} of new covariate values to
    predict with. If not supplied, the original data are used.}
  \item{backTransform}{Logical. If \code{TRUE}, returned values are on the 
    original (e.g. probability, abundance) scale.}
  \item{na.rm}{Logical. Should missing values be removed when predicting from
    original data?}
  \item{appendData}{Logical. Should covariate data used for prediction be appended 
    to the output data frame?}
  \item{level}{The confidence interval to calculate. For example \code{0.95}
    results in lower and upper bounds for a 95\% confidence interval. If set
    to \code{NULL}, no SE or confidence intervals will be returned.}
  \item{re.form}{For \code{unmarkedFit} types that support random effects,
    should the random effects be included in the prediction? If \code{NULL},
    they will be, if \code{NA}, they will not be.}
  \item{se.fit}{Logical. For predicting from \code{unmarkedFitOccuMulti}
    and \code{occuMS} only. Should standard errors and confidence intervals 
    be bootstrapped?}
  \item{species}{Which species (integer or species name as a string) should
    predictions be calculated for? For code{unmarkedFitOccuMulti} only. If
    multiple species are provided, the co-occurence probability is returned.
    See \code{\link{occuMulti}}.}
  \item{cond}{Which species (integer or species name as a string) should
    predictions be calculated conditional on? If conditional on species presence,
    supply just the species name; if species absence, put a minus sign in front
    of the species name, e.g. "-coyote". See \code{\link{occuMulti}}.}
  \item{nsims}{Number of bootstrap simulations to use. Relevant only for some
    predict types.}
  \item{func}{A function to apply to bootstrapped \code{unmarkedRanef} samples;
    see details.}
  \item{...}{Other arguments, currently ignored} 
}

\value{For most methods, a \code{data.frame} with four columns: the
 predicted values, standard errors, and lower and upper bounds. If
 \code{appendData = TRUE} covariate data are also in the output. For
 \code{unmarkedRanef}, an array is returned. See details.
}

\details{
The \code{predict} method for \code{unmarkedRanef} objects generates
the empirical Bayes posterior predictive distribution for functions of the 
random variables (latent abundance or occurrence).

In addition to the output object from \code{ranef}, you must also supply a 
custom function to argument \code{func}. The function must take as input a matrix
with dimensions M x T, where M is the number of sites and T is the number of 
primary periods (T=1 for single-season models). The output of this function should 
be a vector or matrix containing the derived parameters of interest. 

You may also manually set the number of draws from the posterior predictive
distribution with argument \code{nsims}; the default is 100.

The output of \code{predict} will be a vector or array with one more dimension
than the output of the function supplied \code{func}, corresponding to the number
of draws requested \code{nsims}. For example, if \code{func}
outputs a scalar, the output of \code{predict} will be a vector with length
equal to \code{nsims}. If \code{func} outputs a 3x2 matrix, the output of
\code{predict} will be an array with dimensions 3x2x\code{nsims}.
See \code{\link{ranef}} for an example.

Alternatively, you can use the \code{\link{posteriorSamples}} function on the
\code{ranef} output object to obtain the full posterior predictive distribution. 
This is useful if you are having trouble designing your custom function or if 
you want to obtain multiple different derived parameters from the same posterior
predictive distribution.
}
