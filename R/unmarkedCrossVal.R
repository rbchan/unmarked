
setGeneric("crossVal", function(object,
    method=c("Kfold","holdout","leaveOneOut"), folds=10, holdoutPct=0.25,
    statistic=RMSE_MAE, ...) standardGeneric("crossVal"))

setClass("unmarkedCrossVal",
    representation(stats = "data.frame",
                   summary = "data.frame",
                   method = "character",
                   folds = "numeric",
                   holdoutPct = "numeric"),
    validity=function(object){
      errors <- character(0)
      hp <- object@holdoutPct
      if(hp<0|hp>1){
        errors <- c(errors,"holdoutPct must be between 0 and 1")
      }
    }
)

#Constructor of crossVal objects
setMethod("crossVal", "unmarkedFit",
          function(object, method=c("Kfold","holdout","leaveOneOut"),
                   folds=10, holdoutPct=0.25,
                   statistic=RMSE_MAE, parallel=FALSE, ncores, ...){

  method <- match.arg(method, c('Kfold','holdout','leaveOneOut'))

  if(method=="Kfold" & !is.integer(folds) & folds < 0){
    stop("folds must be a positive integer")
  }
  if(method=="holdout" & (holdoutPct>1 | holdoutPct<0)){
    stop("holdoutPct must be a proportion between 0 and 1")
  }

  partitions <- switch(method,
    Kfold = partitionKfold(object, folds=folds),
    holdout = partitionHoldout(object, holdoutPct=holdoutPct),
    leaveOneOut = partitionLeaveOneOut(object)
  )

  n_reps <- length(partitions)

  check_stat <- statistic(object, ...)
  if(!is.numeric(check_stat)||is.null(names(check_stat))){
    stop("Function provided to statistic argument must return a named numeric vector")
  }

  do_crossval <- function(i, object, partitions, statistic, ...){
    newfit <- unmarked::update(object, data=partitions[[i]]$trainData)
    newfit@data <- partitions[[i]]$testData
    if(!is.null(attributes(newfit)$knownOcc)){
      newfit@knownOcc <- rep(FALSE,numSites(newfit@data))
    }

    statistic(newfit, ...)
  }

  if(parallel){
    if(missing(ncores)) ncores <- parallel::detectCores()-1
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    stat_raw <- lapply2(1:n_reps, do_crossval, object,
                                     partitions, statistic, ..., cl = cl)
  } else {
    stat_raw <- lapply2(1:n_reps, do_crossval, object,
                       partitions, statistic, ...)
  }

  stats <- as.data.frame(do.call("rbind", stat_raw))

  summary <- data.frame(Estimate=sapply(stats, mean, na.rm=TRUE),
                        SD=sapply(stats, sd, na.rm=TRUE))

  out <- new("unmarkedCrossVal", stats=stats, summary=summary, method=method,
             folds=folds, holdoutPct=holdoutPct)

  out
})

#Kfold partition function
partitionKfold <- function(object, folds){

  site_inds <- 1:numSites(object@data)
  shuf_site_inds <- sample(site_inds,numSites(object@data))
  fold_inds <- cut(site_inds, breaks=folds, labels=FALSE)

  fold_list <- vector(length=folds,"list")
  for (i in 1:folds){

    trainInds <- shuf_site_inds[fold_inds!=i]
    testInds <- shuf_site_inds[fold_inds==i]

    fold_list[[i]]$trainData <- object@data[trainInds,]
    fold_list[[i]]$testData <-  object@data[testInds,]
  }
  fold_list
}

#Holdout partition function
partitionHoldout <- function(object, holdoutPct){

  site_inds <- 1:numSites(object@data)
  shuf_site_inds <- sample(site_inds,numSites(object@data))

  splitInd <- round(numSites(object@data)*(1-holdoutPct))
  trainInds <- shuf_site_inds[1:splitInd]
  testInds <- shuf_site_inds[(splitInd+1):length(shuf_site_inds)]

  fold_list <- vector(length=1,"list")
  fold_list[[1]]$trainData <- object@data[trainInds,]
  fold_list[[1]]$testData <- object@data[testInds,]

  fold_list
}

#leave-one-out
partitionLeaveOneOut <- function(object){

  fold_list <- vector(length=numSites(object@data),"list")
  for (i in seq_along(fold_list)){
    fold_list[[i]]$trainData <- object@data[-i,]
    fold_list[[i]]$testData <- object@data[i,]
  }
  fold_list

}

setMethod("show", "unmarkedCrossVal", function(object)
{
  st <- object@stats

  if(object@method=='Kfold'){
    cat(paste('Method: k-fold (',object@folds,' folds)\n\n',sep=''))
  } else if(object@method=='holdout'){
    cat(paste('Method: holdout (',round(object@holdoutPct*100),
              '% in test set)\n\n',sep=''))
  } else if(object@method=='leaveOneOut'){
    cat('Method: leave-one-out\n\n')
  }

  for (i in 1:length(st)){
    cat(paste0(names(st)[i],':\n'))
    print(data.frame(object@summary[i,]), row.names=FALSE, digits=4)
    if(i != length(st)) cat('\n')
  }
})

setClass("unmarkedCrossValList",
    representation(stats_list="list",
                   method = "character",
                   folds="numeric",
                   holdoutPct="numeric",
                   sort="character")
)

#CrossVal list constructor
setMethod("crossVal", "unmarkedFitList",
          function(object, method=c("Kfold","holdout","leaveOneOut"),
                   folds=10, holdoutPct=0.25,
                   statistic=RMSE_MAE, parallel=FALSE, ncores,
                   sort = c("none", "increasing", "decreasing"), ...){

    method <- match.arg(method, c('Kfold','holdout','leaveOneOut'))
    sort <- match.arg(sort, c('none','increasing','decreasing'))

    if(missing(ncores)) ncores <- parallel::detectCores()-1
    stats <- lapply(object@fits, crossVal, method, folds,
                    holdoutPct, statistic, parallel, ncores, ...)

    out <- new("unmarkedCrossValList", stats_list=stats, method=method,
               folds=folds, holdoutPct=holdoutPct, sort=sort)

})


setMethod("show", "unmarkedCrossValList", function(object){

  sl <- object@stats_list
  mod_names <- names(sl)
  nfits <- length(sl)
  nstats <- length(sl[[1]]@stats)
  stat_names <- names(sl[[1]]@stats)

  if(object@method=='Kfold'){
    cat(paste('Method: k-fold (',object@folds,' folds)\n\n',sep=''))
  } else if(object@method=='holdout'){
    cat(paste('Method: holdout (',round(object@holdoutPct*100),'% in test set)\n\n',sep=''))
  } else if(object@method=='leaveOneOut'){
    cat('Method: leave-one-out\n\n')
  }

  for (i in 1:nstats){
    cat(paste0(stat_names[i],':\n'))

    stat_sum = lapply(sl, function(x) x@summary[i,])
    stat_sum = do.call("rbind", stat_sum)

    sort_ind <- switch(object@sort,
                       none = 1:nrow(stat_sum),
                       increasing = order(stat_sum$Estimate),
                       decreasing = order(stat_sum$Estimate, decreasing=TRUE))
    stat_sum <- stat_sum[sort_ind, ]

    print(stat_sum, digits=4)
    if(i != nstats) cat('\n')
  }
})

#Function to calculate RMSE and MAE
#Default function for statistic argument
#Returns a named list
RMSE_MAE <- function(object){

  res <- residuals(object)
  if(is.list(res)) res <- unlist(res)

  mae <- mean(abs(res), na.rm=T)
  rmse <- sqrt(mean(res^2, na.rm=T))

  c(`Root mean square error`=rmse, `Mean absolute error`=mae)
}
