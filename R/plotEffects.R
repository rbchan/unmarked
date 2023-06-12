#type_to_covs <- function(umf, type){
#  if(type %in% c("state","psi","lam","lambda","sigma","dist")){
#    return(methods::slot(umf, "siteCovs"))
#  } else if(type %in% c("det","rem","fp","b")){
#    return(methods::slot(umf, "obsCovs"))
#  } else if(type %in% c("phi","transition","col","ext","gamma","omega","iota")){
#    return(methods::slot(umf, "yearlySiteCovs"))
#  }
#  return(NULL)
#}

get_base_newdata <- function(object, type){
  covs <- get_orig_data(object, type)
  out <- lapply(covs, function(x){
    if(is.numeric(x)){
      return(median(x, na.rm=TRUE))
    } else if(is.factor(x)){
      return(factor(levels(x)[1], levels=levels(x)))
    }  else {
      stop("Unknown column type")
    }
  })
  as.data.frame(out)
}

get_cov_seq <- function(covariate, object, type){
  cov_values <- get_orig_data(object, type)[[covariate]]
  if(is.numeric(cov_values)){
    rng <- range(cov_values, na.rm=TRUE)
    return(seq(rng[1], rng[2], length.out=100))
  } else if(is.factor(cov_values)){
    return(factor(levels(cov_values), levels=levels(cov_values)))
  } else {
    stop("Unknown covariate type")
  }
}

setGeneric("plotEffectsData", function(object, ...) standardGeneric("plotEffectsData"))

setMethod("plotEffectsData", "unmarkedFit",
  function(object, type, covariate, level=0.95, ...){

  #umf <- umf_to_factor(object@data)
  nd <- get_base_newdata(object, type)
  if(! covariate %in% names(nd)){
    stop("Covariate not in this submodel", call.=FALSE)
  }
  values <- get_cov_seq(covariate, object, type)
  nd <- nd[rep(1, length(values)),,drop=FALSE]
  nd[[covariate]] <- values

  pr <- predict(object, type=type, newdata=nd, level=level, ...)
  pr$covariate <- covariate
  pr$covariateValue <- values
  pr
})

setGeneric("plotEffects", function(object, ...) standardGeneric("plotEffects"))

setMethod("plotEffects", "unmarkedFit",
  function(object, type, covariate, level=0.95, ...){

  # Get data for plot
  plot_data <- plotEffectsData(object, type, covariate, level, ...)

  # Is covariate a factor?
  is_factor <- is.factor(plot_data$covariateValue)

  # If not a factor (i.e., numeric)
  if(!is_factor){
    # Setup basic plot structure
    plot(plot_data$covariateValue, plot_data$Predicted, type='l',
         xlab=covariate, ylab=paste("Predicted", tolower(object[type]@name)),
         ylim=c(min(plot_data$lower, na.rm=T), max(plot_data$upper, na.rm=T)))
    #Draw error ribbon
    polygon(x=c(plot_data$covariateValue, rev(plot_data$covariateValue)),
            y=c(plot_data$upper, rev(plot_data$lower)),
            border=NA, col='gray90')
    #Draw line on top of ribbon
    lines(plot_data$covariateValue, plot_data$Predicted)

  # If covariate is a factor
  } else {
    # Get number of unique factor levels (groups)
    ngroup <- nrow(plot_data)
    # Calculate width of error bar tips based on # of factor levels
    bw <- ngroup/10

    # Add basic plot of points
    plot(1:ngroup, plot_data$Predicted,
         xlab=covariate, ylab=paste("Predicted", tolower(object[type]@name)),
         xlim=c(1-bw/2,ngroup+bw/2),
         xaxt='n',
         ylim=c(min(plot_data$lower, na.rm=T), max(plot_data$upper, na.rm=T)))

    # Draw new x-axis with factor labels
    axis(1, at=1:ngroup, labels=plot_data$covariateValue)

    # Draw error bars
    for (i in 1:ngroup){
      segments(i, plot_data$lower[i], i, plot_data$upper[i])
      segments(i-bw/2, plot_data$lower[i], i+bw/2, plot_data$lower[i])
      segments(i-bw/2, plot_data$upper[i], i+bw/2, plot_data$upper[i])
    }
    #Draw points on top of bars
    points(1:ngroup, plot_data$Predicted, pch=19)
  }

})
