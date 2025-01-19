setClass("unmarkedFrameOccuRNMulti",
         contains="unmarkedFrameOccuMulti")

setClass("unmarkedFitOccuRNMulti",
         representation(
            detformulas = "list",
            stateformulas = "list",
            modelOccupancy = "numeric"),
         contains = "unmarkedFit")

occuRNMulti <- function(detformulas, stateformulas, data, modelOccupancy,
                        K = rep(10, length(stateformulas)),
                        starts, method="BFGS", se = TRUE, threads=1, ...){
 
  S <- length(stateformulas)
  dep <- create_dep_matrix(stateformulas)

  # Reorder everything to match
  stopifnot(all(names(data@ylist) %in% names(stateformulas)))
  stopifnot(all(names(data@ylist) %in% names(detformulas)))
  species_order <- colnames(dep)
  stateformulas <- stateformulas[species_order]
  stateformulas <- lapply(stateformulas, function(x){
    if(is.call(x)) return(x)
    x[order(match(names(x), c("", species_order)))]
  })
  detformulas <- detformulas[species_order]
  data@ylist <- data@ylist[species_order]

  modocc <- rep(0, S)
  names(modocc) <- names(data@ylist)
  if(!missing(modelOccupancy)){
    stopifnot(all(modelOccupancy %in% names(modocc)))
    modocc[modelOccupancy] <- 1
    not_allowed_occ <- colSums(dep)
    sp_mismatch <- modocc & not_allowed_occ
    if(any(sp_mismatch)){
      stop(paste("Species", 
                 paste(names(modocc)[sp_mismatch], collapse=", "),
                 "cannot use occupancy model."), call.=FALSE)
    }
  }

  data <- as(data, "unmarkedFrameOccuRNMulti")
  gd <- getDesign(data, detformulas, stateformulas)

  if(length(K) == 1){
    K <- rep(K, S)
  }
  stopifnot(length(K) == S)
  Kmin <- sapply(gd$ylist, function(x){
                   apply(x, 1, function(y){
                           if(all(is.na(y))) return(0)
                           max(y, na.rm=TRUE)
                    })
                })

  nP <- max(gd$det_ind)
  
  has_occ <- sum(modocc) > 0
  if(has_occ){
    sp_occ <- names(modocc)[modocc == 1]
    sp_abun <- names(modocc)[modocc == 0]
    occ_ind <- gd$state_ind[grepl(paste(paste0("^",sp_occ), collapse="|"),
                                  rownames(gd$state_ind)),,drop=FALSE]
    occ_rng <- min(occ_ind):max(occ_ind)
    abun_ind <- gd$state_ind[grepl(paste(paste0("^",sp_abun), collapse="|"),
                                  rownames(gd$state_ind)),,drop=FALSE]
    abun_rng <- min(abun_ind):max(abun_ind)
  } else {
    state_rng <- min(gd$state_ind):max(gd$state_ind)
  }

  det_rng <- min(gd$det_ind):max(gd$det_ind)
  
  if(missing(starts)) starts <- rep(0, nP)
  names(starts) <- gd$par_names

  nll_C <- function(pars){
    nll_occuRNMulti(pars, gd$state_ind-1, gd$det_ind-1, S, modocc,
                    gd$ylist, gd$state_dm, gd$det_dm, dep, K, Kmin,
                    gd$miss, gd$site_miss, threads)
  }

  fm <- optim(starts, nll_C, method=method, hessian=se, ...)
  covMat <- invertHessian(fm, nP, se)

  fmAIC <- 2 * fm$value + 2 * nP
  #----------------------------------------------------------------------------

  #Format output---------------------------------------------------------------
  ests <- fm$par
  names(ests) <- gd$par_names

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates = ests[det_rng],
                          covMat = as.matrix(covMat[det_rng, det_rng]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  if(has_occ){
    lam <- unmarkedEstimate(name = "Abundance", short.name = "lam",
                            estimates = ests[abun_rng],
                            covMat = as.matrix(covMat[abun_rng, abun_rng]),
                            invlink = "exp",
                            invlinkGrad = "exp")
    psi <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                          estimates = ests[occ_rng],
                          covMat = as.matrix(covMat[occ_rng, occ_rng]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

    estimateList <- unmarkedEstimateList(list(lam=lam, psi=psi, det=det))
  } else {
    lam <- unmarkedEstimate(name = "Abundance", short.name = "lam",
                            estimates = ests[state_rng],
                            covMat = as.matrix(covMat[state_rng, state_rng]),
                            invlink = "exp",
                            invlinkGrad = "exp")
    estimateList <- unmarkedEstimateList(list(lam=lam, det=det))
  }

  umfit <- new("unmarkedFitOccuRNMulti", fitType = "occuRNMulti", call = match.call(),
                detformulas = detformulas, stateformulas = stateformulas,
                modelOccupancy = modocc,
                formula = ~1, data = data,
                #sitesRemoved = designMats$removed.sites,
                estimates = estimateList, AIC = fmAIC, opt = fm,
                negLogLike = fm$value, nllFun = nll_C)

  return(umfit)
}

setMethod("getDesign", "unmarkedFrameOccuRNMulti",
  function(umf, detformulas, stateformulas){
  M <- nrow(umf@ylist[[1]])
  J <- ncol(umf@ylist[[1]])
  
  ylist <- umf@ylist

  # Site covs
  sc <- umf@siteCovs
  if(is.null(sc)) sc <- data.frame(dummy_=rep(1, M))

  mm <- lapply(stateformulas, function(s){
    if(!is.list(s)){
      mf <- model.frame(s, sc, na.action=stats::na.pass)
      return(model.matrix(s, mf))
    }
    lapply(s, function(x){
           mf <- model.frame(x, sc, na.action=stats::na.pass)
           model.matrix(x, mf)
    })
  })
  all_mm <- lapply(mm, function(x){
                     if(is.matrix(x)) return(x)
                     do.call(cbind, x)
             })
  all_mm <- do.call(cbind, all_mm)
  miss_site_covs <- apply(all_mm, 1, function(x) any(is.na(x)))

  # Obs covs
  sc_long_ind <- rep(1:M, each=J)
  sc_long <- sc[sc_long_ind,]
  if(is.null(umf@obsCovs)){
    oc <- sc_long
  } else {
    oc <- cbind(sc_long, umf@obsCovs)
  }

  vv <- lapply(detformulas, function(x){
    mf <- model.frame(x, oc, na.action=stats::na.pass)
    model.matrix(x, mf)
  })
  all_vv <- do.call(cbind, vv)
  miss_obs_covs <- apply(all_vv, 1, function(x) any(is.na(x)))
  miss_obs_covs <- matrix(miss_obs_covs, M, J, byrow=TRUE)

  # Indices
  state_n <- unlist(lapply(mm, function(x){
    if(is.list(x)){
      lapply(x, function(y) ncol(y))
    } else {
      ncol(x)
    }
  }))
  state_end <- cumsum(state_n)
  state_start <- state_end - state_n + 1
  state_ind <- cbind(state_start, state_end)

  state_prefix <- rep(rownames(state_ind), state_n)
  state_prefix <- gsub(".", ":", state_prefix, fixed=TRUE)
  state_prefix <- paste0("[",state_prefix,"]")
  
  state_par <- unlist(lapply(mm, function(x){
    if(is.matrix(x)) return(colnames(x))
    lapply(x, function(y){
      return(colnames(y))
    })
  }))
  stopifnot(length(state_prefix) == length(state_par))
  state_par <- paste(state_prefix, state_par)
  
  det_n <- sapply(vv, ncol)
  det_end <- cumsum(det_n) + max(state_end)
  det_start <- det_end - det_n + 1
  det_ind <- cbind(det_start, det_end)

  det_rng <- min(det_ind):max(det_ind)

  det_prefix <- paste0("[",rep(rownames(det_ind), det_n),"]")
  det_par <- unlist(lapply(vv, colnames))
  stopifnot(length(det_prefix) == length(det_par))
  det_par <- paste(det_prefix, det_par)
  par_names <- c(state_par, det_par)

  # Missing values
  miss <- (is.na(umf@ylist[[1]]) | miss_obs_covs) * 1
  site_miss <- (apply(umf@ylist[[1]], 1, function(x) all(is.na(x))) | miss_site_covs) * 1

  list(par_names = par_names, state_ind = state_ind, det_ind = det_ind,
              ylist = ylist, state_dm = mm, det_dm = vv,
              miss = miss, site_miss = site_miss)
})

create_dep_matrix <- function(stateformulas){
  S <- length(stateformulas)
  dep <- matrix(0, S, S)
  rownames(dep) <- colnames(dep) <- names(stateformulas)

  # Dominant species
  is_dom <- sapply(stateformulas, is.call)
  if(sum(is_dom) == 0){
    stop("Must be at least one dominant species not depending on any others",
         call.=FALSE)
  }
  allowed_species <- sort(names(stateformulas)[is_dom])

  # 2nd level
  is_2nd <- sapply(stateformulas, function(x){
    if(is.call(x)) return(FALSE)
    stopifnot(is.list(x))
    nm <- names(x)[names(x) != ""]
    all(nm %in% allowed_species)
  })
  species_2nd <- sort(names(stateformulas)[is_2nd])
  allowed_species <- c(allowed_species, species_2nd)
  
  for (i in species_2nd){
    dep_sp <- names(stateformulas[[i]])
    dep_sp <- dep_sp[dep_sp != ""]
    dep[i, dep_sp] <- 1
  }

  # 3rd level
  is_3rd <- sapply(1:length(stateformulas), function(i){
    if(names(stateformulas)[i] %in% allowed_species) return(FALSE)
    x <- stateformulas[[i]]
    stopifnot(is.list(x))
    nm <- names(x)[names(x) != ""]
    all(nm %in% allowed_species)
  })
  names(is_3rd) <- names(stateformulas)
  
  species_3rd <- sort(names(stateformulas)[is_3rd])
  allowed_species <- c(allowed_species, species_3rd)
  
  for (i in species_3rd){
    dep_sp <- names(stateformulas[[i]])
    dep_sp <- dep_sp[dep_sp != ""]
    dep[i, dep_sp] <- 1
  }

  if(any(!names(stateformulas) %in% allowed_species)){
    stop("Invalid stateformula list", call.=FALSE)
  }
  
  dep <- dep[allowed_species, allowed_species]

  dep
}

# Compute linear combinations of estimates in unmarkedFit objects.
setMethod("linearComb",
    signature(obj = "unmarkedFitOccuRNMulti", coefficients = "matrixOrVector"),
    function(obj, coefficients, type, offset = NULL, re.form=NULL)
{
  stop("This method not supported for occuRNMulti fits, use predict instead.",
       call.=FALSE)
})

setMethod("predict", "unmarkedFitOccuRNMulti",
          function(object, type, species, newdata,
                   level = 0.95, nsims = 100, ...){

  if(!missing(species)){
    stopifnot(species %in% names(object@data@ylist))
  }

  se <- TRUE
  if(is.null(hessian(object))){
    se = FALSE
  }

  if(type == "state"){

    depmat <- unmarked:::create_dep_matrix(object@stateformulas)

    top <- apply(depmat, 1, function(x) sum(x) == 0)
    sp_top <- colnames(depmat)[top]
    lev2 <- apply(depmat[!top,!top,drop=FALSE], 1, function(x) sum(x) == 0)
    sp_2nd <- names(lev2)[lev2]
    sp_3rd <- names(lev2)[!lev2]

    if(missing(newdata) || is.null(newdata)){
      newdata <- NULL 
      X <- lapply(object@stateformulas, function(x){
        if(is.list(x)){
          lapply(x, function(z){
                   mf <- model.frame(z, object@data@siteCovs, na.action=stats::na.pass)
                   model.matrix(z, mf)
          })
        } else {
          mf <- model.frame(x, object@data@siteCovs, na.action=stats::na.pass)
          model.matrix(x, mf)
        }
      })
      nr <- nrow(object@data@siteCovs)
    } else {
      X <- lapply(object@stateformulas, function(x){
        if(is.list(x)){
          lapply(x, function(z) unmarked:::make_mod_matrix(z, object@data@siteCovs, newdata)$X)
        } else {
          unmarked:::make_mod_matrix(x, object@data@siteCovs, newdata)$X
        }
      })
      nr <- nrow(newdata)
    }
    
    coef_est <- matrix(coef(object), nrow=1)
    colnames(coef_est) <- names(coef(object))
    point_est <- calc_dependent_response(coef_est, X, object, nr, sp_top, sp_2nd, sp_3rd, newdata)


    if(se & !is.null(level)){
      message('Bootstrapping confidence intervals with ',nsims,' samples')
      samps <- MASS::mvrnorm(nsims, mu=coef(object), Sigma = vcov(object))
      post <- calc_dependent_response(samps, X, object, nr, sp_top, sp_2nd, sp_3rd, newdata)

      # Summarize bootstrap
      cis <- lapply(post, function(x){
        data.frame(
          SE = apply(x, 1, sd, na.rm=TRUE),
          lower = apply(x, 1, quantile, (1-level)/2, na.rm=TRUE),
          upper = apply(x, 1, quantile, 1 - (1-level)/2, na.rm=TRUE)
        )
      })
    } else {
      cis <- lapply(point_est, function(x){
        data.frame(SE = rep(NA, nrow(x)), lower=NA, upper=NA)
      })
    }

    out <- mapply(function(x, y){
      cbind(Predicted = x, y)
    }, point_est, cis, SIMPLIFY=FALSE)

    if(!missing(species)){
      out <- out[[species]]
    }
    return(out)
  } else if(type == "det"){

    det_species <- names(object@data@ylist)
    if(!missing(species)){
      det_species <- species
    }

    gd <- unmarked:::getDesign(object@data, object@detformulas, object@stateformulas)
    od <- unmarked:::get_orig_data(object, "det")
    chunk_size <- 70
  
    if(missing(newdata)) newdata <- NULL
    out <- lapply(det_species, function(i){
      form <- object@detformulas[[i]]
      cf <- coef(object)
      inds <- which(grepl(paste0("p([",i,"]"), names(cf), fixed=TRUE))
      new_est <- object@estimates@estimates$det
      new_est@estimates <- cf[inds]
      new_est@fixed <- 1:length(inds)
      if(se){
        new_est@covMat <- vcov(object)[inds,inds,drop=FALSE]
        new_est@covMatBS <- object@covMatBS[inds,inds,drop=FALSE]
      } else {
        new_est@covMat <- matrix(NA, nrow=length(inds), ncol=length(inds))
        new_est@covMatBS <- matrix(NA, nrow=length(inds), ncol=length(inds))
      }
     
      if(!is.null(newdata)){
        X <- unmarked:::make_mod_matrix(form, od, newdata, re.form=NULL)$X
      } else {
        X <- gd$det_dm[[i]] 
      }
      nr <- nrow(X)
      sep <- rep(1:ceiling(nr/chunk_size), each=chunk_size, length.out=nr)

      x_chunk <- lapply(unique(sep),
                    function(i) as.matrix(X[sep==i,,drop=FALSE]))

      prmat <- lapply(x_chunk, function(x_i){
        has_na <- apply(x_i, 1, function(x_i) any(is.na(x_i)))
        # Work around linearComb bug where there can't be NAs in inputs
        x_i[has_na,] <- 0
        lc <- linearComb(new_est, x_i)
        lc <- backTransform(lc)
        out <- data.frame(Predicted=coef(lc), SE=NA, lower=NA, upper=NA)
        if(!is.null(level)){
          se <- SE(lc)
          ci <- confint(lc, level=level)
          out$SE <- se
          out$lower <- ci[,1]
          out$upper <- ci[,2]
        }
        out[has_na,] <- NA
        out
      })
      prmat <- do.call(rbind, prmat)
      rownames(prmat) <- NULL
      prmat
    })
    names(out) <- det_species
    if(length(out) == 1) out <- out[[1]]
    return(out)
  } else {
    stop("type must be state or det", call.=FALSE)
  }
})

# Function to extract correct parts of state coefs vector
get_species_b <- function(b, object, sp){
  b <- b[!grepl("^p\\(", names(b))] # remove detection coefficients
  nm <- names(b)
  allsp <- names(object@data@ylist)
  ind <- which(allsp == sp)
  out <- list()
  other <- c(1:length(allsp))[-ind]
  match1 <- grepl(paste0("[", allsp[ind], "]"), nm, fixed=TRUE)
  out[[1]] <- b[match1]

  match2 <- grepl(paste0("[", sp, ":", allsp[other[1]],"]"), nm, fixed=TRUE)
  if(sum(match2 > 0)){
    out[[allsp[other[1]]]] <- b[match2]
  }
  match3 <- grepl(paste0("[", sp, ":", allsp[other[2]],"]"), nm, fixed=TRUE)
  if(sum(match3 > 0)){
    out[[allsp[other[2]]]] <- b[match3]
  }
  out
}

calc_dependent_response <- function(samps, X, object, nr, sp_top, sp_2nd, sp_3rd, 
                                    newdata){
  nsims <- nrow(samps)

  all_species <- names(object@data@ylist)

  post <- lapply(1:length(all_species), function(x) matrix(NA, nr, nsims))
  names(post) <- all_species

  for (k in 1:nsims){

    b <- samps[k,]
    ests <- list()

    # Top-level species that don't depend on other species abundance
    for (i in sp_top){
      beta <- get_species_b(b, object, i)[[1]]
      ests[[i]] <- exp(X[[i]] %*% beta)
      post[[i]][,k] <- ests[[i]]
    }
    # 2nd-level species that depend only on top-level species
    for (i in sp_2nd){
      if(object@modelOccupancy[i]){
        ilink <- plogis
      } else {
        ilink <- eval(str2lang(object@estimates@estimates$lam@invlink))
      }
      Xsub <- X[[i]]
      beta <- get_species_b(b, object, i)
      lp <- Xsub[[1]] %*% beta[[1]]
      for (j in 2:length(Xsub)){
        other_sp <- names(Xsub)[j]

        # Use fixed species abundances if they are in newdata
        if(other_sp %in% names(newdata)){
          other_sp_abun <- newdata[[other_sp]]
        } else {
          other_sp_abun <- ests[[other_sp]]
        }

        add <- Xsub[[other_sp]] %*% beta[[other_sp]] * other_sp_abun
        lp <- lp + add
      }
      ests[[i]] <- ilink(lp)
      post[[i]][,k] <- ests[[i]]
    }
    # 3rd-level species that depend on at least one 2nd-level species
    for (i in sp_3rd){
      if(object@modelOccupancy[i]){
        ilink <- plogis
      } else {
        ilink <- eval(str2lang(object@estimates@estimates$lam@invlink))
      }
      Xsub <- X[[i]]
      beta <- get_species_b(b, object, i)
      lp <- Xsub[[1]] %*% beta[[1]]
      for (j in 2:length(Xsub)){
        other_sp <- names(Xsub)[j]

        # Use fixed species abundances if they are in newdata
        if(other_sp %in% names(newdata)){
          other_sp_abun <- newdata[[other_sp]]
        } else {
          other_sp_abun <- ests[[other_sp]]
        }

        add <- Xsub[[other_sp]] %*% beta[[other_sp]] * other_sp_abun
        lp <- lp + add
      }
      ests[[i]] <- ilink(lp)
      post[[i]][,k] <- ests[[i]]
    }
  }
  post
}


setMethod("getP", "unmarkedFitOccuRNMulti", function(object)
{

  ylist <- object@data@ylist
  S <- length(ylist)
  M <- nrow(ylist[[1]])
  J <- ncol(ylist[[1]])

  pr <- predict(object, type = "det", level=NULL)
  pr <- lapply(pr, function(x){
    matrix(x$Predicted, nrow=M, ncol=J, byrow=TRUE)
  })
  pr
})

setMethod("fitted", "unmarkedFitOccuRNMulti", function(object, na.rm=FALSE){
  S <- length(object@data@ylist)
  M <- nrow(object@data@ylist[[1]])
  J <- ncol(object@data@ylist[[1]])

  # This will be lambda for RN models and psi for occupancy models
  state <- predict(object, type="state", level=NULL)
  # Reorganize into M x J matrices
  state <- lapply(state, function(x){
    matrix(rep(x$Predicted, J), nrow=M, ncol=J)
  })
  # This will be r for RN models and p for occupancy models
  p <- getP(object)
  
  # Need to calculate this differently depending on if state is abundance or occupancy
  fitted <- lapply(1:S, function(i){
    if(object@modelOccupancy[i]){
      # If occupancy, fitted is psi * p
      state[[i]] * p[[i]]
    } else {
      # If abundance
      1 - exp(-state[[i]]*p[[i]]) ## analytical integration.
    }
  })
  names(fitted) <- names(object@data@ylist)
  fitted
})

setMethod("residuals", "unmarkedFitOccuRNMulti", function(object, ...){
  ylist <- object@data@ylist
  S <- length(ylist)
  ft <- fitted(object)
  res <- lapply(1:S, function(i){
    ylist[[i]] - ft[[i]] 
  })
  names(res) <- names(ylist)
  res
})

setMethod("SSE", "unmarkedFitOccuRNMulti", function(fit, ...){
    r <- do.call(rbind, residuals(fit))
    return(c(SSE = sum(r^2, na.rm=T)))
})

setMethod("nonparboot", "unmarkedFitOccuRNMulti",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
  stop("Method not supported for this fit type", call.=FALSE)
})

setMethod("simulate", "unmarkedFitOccuRNMulti", 
  function(object, nsim = 1, seed = NULL, na.rm=TRUE){

  ylist <- object@data@ylist
  M <- nrow(ylist[[1]])
  J <- ncol(ylist[[1]])
  S <- length(ylist)
  occ <- object@modelOccupancy

  state <- predict(object, type="state", level=NULL)
  state <- lapply(state, function(x) x$Predicted)
  p <- getP(object)

  sims <- lapply(1:nsim, function(i){
    out <- lapply(1:S, function(s){
      y <- matrix(NA, M, J)
      if(occ[s]){
        z <- rbinom(M, 1, state[[s]])
        for (j in 1:J){
          y[,j] <- rbinom(M, 1, p[[s]][,j])
        }
        y
      } else {
        N <- rpois(M, state[[s]])
        for (j in 1:J){
          y[,j] <- rbinom(M, N, p[[s]][,j])
        }
        y[y>0] <- 1
        y
      }
    })
    names(out) <- names(ylist)
    out
  })
  
  sims
})

setMethod("replaceY", "unmarkedFrameOccuRNMulti",
          function(object, newY, replNA=TRUE, ...){
      if(replNA){
        newY <- mapply(function(x, y){ is.na(x) <- is.na(y); x},
                       newY , object@ylist, SIMPLIFY=FALSE)
      }
      object@ylist <- newY
      object
})


setMethod("update", "unmarkedFitOccuRNMulti",
    function(object, detformulas, stateformulas, ..., evaluate = TRUE)
{

    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    if(!missing(detformulas)){
      call[["detformulas"]] <- detformulas
    } else {
      call[["detformulas"]] <- object@detformulas
    }
    if(!missing(stateformulas)){
      call[["stateformulas"]] <- stateformulas
    } else {
      call[["stateformulas"]] <- object@stateformulas
    }
    extras <- match.call(call=sys.call(-1),
                         expand.dots = FALSE)$...
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing])
            call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
            }
        }
    if (evaluate)
        eval(call, parent.frame(2))
    else call
})

setMethod("ranef", "unmarkedFitOccuRNMulti",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
   stop("Not currently supported for occuRNMulti", call.=FALSE)
})
