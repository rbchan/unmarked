
setClassUnion("matrixOrVector", c("matrix","numeric"))

# Class to store actual parameter estimates
setClass("unmarkedEstimate",
    representation(
        name = "character",
		short.name = "character",
        estimates = "numeric",
        covMat = "matrix",
        fixed = "numeric",
        covMatBS = "optionalMatrix",
        invlink = "character",
        invlinkGrad = "character",
        randomVarInfo= "list"),
    validity = function(object) {
        errors <- character(0)
        if(nrow(object@covMat) != length(object@estimates)) {
        errors <- c(errors,
            "Size of covMat does not match length of estimates.")
        }
    if(length(errors) > 0)
        errors
    else
        TRUE
    })

setClass("unmarkedEstimateList",
    representation(estimates = "list"),
    validity = function(object) {
        errors <- character(0)
        for(est in object@estimates) {
            if(!is(est, "unmarkedEstimate")) {
                errors <- c("At least one element of unmarkedEstimateList is not an unmarkedEstimate.")
                break
                }
            }
        if(length(errors) == 0) {
            return(TRUE)
        } else {
            return(errors)
            }
    })

setMethod("show", "unmarkedEstimateList",
    function(object) {
      for(est in object@estimates) {
        show(est)
        cat("\n")
      }
    })

setMethod("summary", "unmarkedEstimateList",
    function(object)
{
    sumList <- list()
    for(i in 1:length(object@estimates)) {
        sumList[[i]] <- summary(object@estimates[[i]])
        cat("\n")
        }
    names(sumList) <- names(object@estimates)
    invisible(sumList)
})

setGeneric("estimates",
    function(object) {
      standardGeneric("estimates")
    })

setMethod("estimates", "unmarkedEstimate",
    function(object) {
      object@estimates
    })

setMethod("estimates", "unmarkedEstimateList",
		function(object) {
			object@estimates
		})

unmarkedEstimateList <- function(l) {
  new("unmarkedEstimateList", estimates = l)
}



unmarkedEstimate <- function(name, short.name, estimates, covMat, fixed=NULL,
                             invlink, invlinkGrad, randomVarInfo=list())
{

    if(is.null(fixed)) fixed <- 1:length(estimates)
    new("unmarkedEstimate",
        name = name,
        short.name = short.name,
        estimates = estimates,
        covMat = covMat,
        fixed = fixed,
        invlink = invlink,
        invlinkGrad = invlinkGrad,
        randomVarInfo = randomVarInfo)
}


setMethod("show", signature(object = "unmarkedEstimate"),
    function(object)
{

    has_random <- methods::.hasSlot(object, "randomVarInfo") &&
      length(object@randomVarInfo) > 0

    fixed <- 1:length(object@estimates)
    if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
    ests <- object@estimates[fixed]
    SEs <- SE(object)[fixed]
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <- !(length(ests) == 1 |
                       identical(names(ests), "(Intercept)") |
                       identical(names(ests), "1"))

    cat(object@name,":\n", sep="")

    if(has_random){
      print_randvar_info(object@randomVarInfo)
      cat("\nFixed effects:\n")
    }

    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)

})



setMethod("summary", signature(object = "unmarkedEstimate"),
    function(object)
{
    fixed <- 1:length(object@estimates)
    if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
    ests <- object@estimates[fixed]
    SEs <- SE(object)[fixed]
    Z <- ests/SEs
    p <- 2*pnorm(abs(Z), lower.tail = FALSE)
    printRowNames <-
        !(length(ests) == 1 | identical(names(ests), "(Intercept)") | identical(names(ests), "1"))
    invlink <- object@invlink
    link <- switch(invlink,
                   exp = "log",
                   logistic = "logit",
                   cloglog = "cloglog")
    cat(object@name, " (", link, "-scale)", ":\n", sep="")

    has_random <- methods::.hasSlot(object, "randomVarInfo") && 
      length(object@randomVarInfo) > 0
    if(has_random){
      print_randvar_info(object@randomVarInfo)
      cat("\nFixed effects:\n")
    }

    outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                        check.names = FALSE)
    print(outDF, row.names = printRowNames, digits = 3)
    invisible(outDF)
})



# Compute linear combinations of estimates in unmarkedEstimate objects.

setMethod("linearComb",
    signature(obj = "unmarkedEstimate", coefficients = "matrixOrVector"),
    function(obj, coefficients, offset = NULL, re.form = NULL)
{
    if(!is(coefficients, "matrix"))
        coefficients <- t(as.matrix(coefficients))
    est <- obj@estimates
    covMat <- obj@covMat
    if(!is.null(re.form) & .hasSlot(obj, "fixed")){
      est <- est[obj@fixed]
      covMat <- covMat[obj@fixed, obj@fixed, drop=FALSE]
    }
    stopifnot(ncol(coefficients) == length(est))
    if (is.null(offset))
        offset <- rep(0, nrow(coefficients))
    e <- as.vector(coefficients %*% est) + offset
    v <- coefficients %*% covMat %*% t(coefficients)
    if (!is.null(obj@covMatBS)) {
        v.bs <- coefficients %*% obj@covMatBS %*% t(coefficients)
    } else {
        v.bs <- NULL
    }
    umelc <- new("unmarkedLinComb", parentEstimate = obj,
                 estimate = e, covMat = v, covMatBS = v.bs,
                 coefficients = coefficients)
    umelc
})


# Transform an unmarkedEstimate object to it's natural scale.
#setMethod("backTransform",
#    signature(obj = "unmarkedEstimate"),
#    function(obj) {
#      stopifnot(length(obj@estimates) == 1)
#      e <- eval(call(obj@invlink,obj@estimates))
#      v <- (eval(call(obj@invlinkGrad,obj@estimates)))^2 * obj@covMat
#
#      if(is(obj, "unmarkedEstimateLinearComb")) {
#        coef <- obj@coefficients
#        orig <- obj@originalEstimate
#      } else {
#        coef <- 1
#        orig <- obj
#      }
#
#      umebt <- new("unmarkedEstimateBackTransformed",
#          name = paste(obj@name,"transformed to native scale"),
#          estimates = e, covMat = v,
#          invlink = "identity", invlinkGrad = "identity",
#          originalEstimate = orig, coefficients = coef,
#          transformation = obj@invlink)
#      umebt
#    })



# backTransform is only valid for an unmarkedEstimate of length = 1.
# can backtranform a fit directly if it has length 1
# o.w. give error
setMethod("backTransform", "unmarkedEstimate", function(obj)
{
    if(length(obj@estimates) == 1) {
        lc <- linearComb(obj, 1)
        return(backTransform(lc))
    } else {
        stop("Cannot directly back-transform an unmarkedEstimate with length > 1.\nUse linearComb() and then backTransform() the resulting scalar linear combination.")
    }
})


# Compute standard error of an unmarkedEstimate object.
setMethod("SE", signature(obj = "unmarkedEstimate"), function(obj, fixedOnly=TRUE)
{
    sqrt(diag(vcov(obj, fixedOnly=fixedOnly)))
})


setMethod("[", signature("unmarkedEstimateList"),
    function(x, i, j, drop)
{
    x@estimates[[i]]
})


setMethod("names", "unmarkedEstimateList",
    function(x)
{
    names(x@estimates)
})


setMethod("coef", "unmarkedEstimate",
    function(object, altNames = TRUE, fixedOnly=TRUE, ...)
{
    coefs <- object@estimates
    if(fixedOnly){
      fixed <- 1:length(coefs)
      if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
      coefs <- coefs[fixed]
    }
    names(coefs)[names(coefs) == "(Intercept)"] <- "Int"
    if(altNames) {
        names(coefs) <- paste(object@short.name, "(", names(coefs), ")",
                              sep="")
    }
    coefs
})


setMethod("vcov", "unmarkedEstimate",
    function(object, fixedOnly=TRUE, ...)
{
        v <- object@covMat
        if(fixedOnly){
          fixed <- 1:nrow(v)
          if(methods::.hasSlot(object, "fixed")) fixed <- object@fixed
          v <- as.matrix(v[fixed,fixed])
        }

        rownames(v) <- colnames(v) <- names(coef(object, fixedOnly=fixedOnly))
        v
})


setMethod("confint", "unmarkedEstimate",
    function(object, parm, level = 0.95)
{
    if(missing(parm)) parm <- object@fixed
    ests <- object@estimates[parm]
    ses <- SE(object)[parm]
    z <- qnorm((1-level)/2, lower.tail = FALSE)
    lower.lim <- ests - z*ses
    upper.lim <- ests + z*ses
    ci <- as.matrix(cbind(lower.lim, upper.lim))
    rownames(ci) <- names(coef(object))[parm]
    colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
    ci
})





