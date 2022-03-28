get_xlev <- function(data, model_frame){
  fac_col <- data[, sapply(data, is.factor), drop=FALSE]
  xlevs <- lapply(fac_col, levels)
  xlevs[names(xlevs) %in% names(model_frame)]
}

get_reTrms <- function(formula, data, newdata=NULL){
  fb <- lme4::findbars(formula)
  mf <- model.frame(lme4::subbars(formula), data, na.action=stats::na.pass)
  if(is.null(newdata)) return(lme4::mkReTrms(fb, mf))
  new_mf <- model.frame(stats::terms(mf), newdata, na.action=stats::na.pass,
                        xlev=get_xlev(data, mf))
  lme4::mkReTrms(fb, new_mf, drop.unused.levels=FALSE)
}

get_Z <- function(formula, data, newdata=NULL){
  if(is.null(lme4::findbars(formula))){
    if(is.null(newdata)){
      return(Matrix::Matrix(matrix(0, nrow=nrow(data), ncol=0),sparse=TRUE))
    } else{
      return(Matrix::Matrix(matrix(0, nrow=nrow(newdata), ncol=0),sparse=TRUE))
    }
  }
  check_formula(formula, data)
  Zt <- get_reTrms(formula, data, newdata)$Zt
  Z <- t(as.matrix(Zt))
  Matrix::Matrix(Z, sparse=TRUE)
}

get_group_vars <- function(formula){
  rand <- lme4::findbars(formula)
  ifelse(is.null(rand), 0, length(rand))
}

get_nrandom <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(as.array(0))

  out <- sapply(rand, function(x){
    col_nm <- as.character(x[[3]])
    length(unique(data[[col_nm]]))
  })
  as.array(out)
}

has_random <- function(formula){
  length(lme4::findbars(formula)) > 0
}

sigma_names <- function(formula, data){
  if(!has_random(formula)) return(NA_character_)
  nms <- get_reTrms(formula, data)$cnms
  nms <- paste0(unlist(nms), "|", names(nms))
  nms <- gsub("(Intercept)", "1", nms, fixed=TRUE)
  #paste0("sigma [", nms, "]")
  nms
}

check_formula <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(is.null(rand)) return(invisible())

  char <- paste(formula, collapse=" ")
  if(grepl(":|/", char)){
    stop("Nested random effects (using / and :) are not supported",
         call.=FALSE)
  }
  theta <- get_reTrms(formula, data)$theta
  if(0 %in% theta){
    stop("Correlated slopes and intercepts are not supported. Use || instead of |.",
         call.=FALSE)
  }
}

split_formula <- function(formula){
  if(length(formula) != 3) stop("Double right-hand side formula required")
  char <- lapply(formula, function(x){
            paste(deparse(x), collapse="")
          })
  p1 <- as.formula(char[[2]])
  p2 <- as.formula(paste("~", char[[3]]))
  list(p1, p2)
}

is_tmb_fit <- function(mod){
  if(!methods::.hasSlot(mod, "TMB")) return(FALSE)
  !is.null(mod@TMB)
}

get_b_vector <- function(tmb_report, type){
  #sdr <- TMB::sdreport(TMB)
  bname <- paste0("b_",type)
  bpar <- tmb_report$par.random
  bpar <- bpar[grepl(bname, names(bpar))]
  if(length(bpar)==0) return(NULL)
  bpar
}

get_joint_cov <- function(tmb_report, type=NULL, remove_sigma=TRUE){
  full <- tmb_report$jointPrecision
  if(is.null(full)){
    full <- tmb_report$cov.fixed
    colnames(full) <- rownames(full) <- names(tmb_report$par)
  } else {
    full <- solve(full)
  }

  if(is.null(type)) return(full)
  keep <- grepl(paste0("_",type), colnames(full))
  out <- full[keep,keep,drop=FALSE]
  if(!remove_sigma) return(out)
  keep2 <- !grepl(paste0("lsigma_",type), colnames(out))
  out <- out[keep2,keep2,drop=FALSE]
}

tmbfit_has_random <- function(mod, type){
  paste0(type,"RE") %in% names(mod@estimates@estimates)
}

use_tmb_bootstrap <- function(mod, type, re.form){
  is.null(re.form) && is_tmb_fit(mod) && tmbfit_has_random(mod, type)
}

# Gather information about grouping variables for a given submodel
get_randvar_info <- function(tmb_report, type, formula, data){
  ngv <- get_group_vars(formula)
  if(ngv == 0) return(list()) #Return blank list if there are no grouping variables

  sigma_type <- paste0("lsigma_",type)
  sigma_ind <- grepl(sigma_type, get_fixed_names(tmb_report))
  sigma_est <- tmb_report$par.fixed[sigma_ind]
  sigma_cov <- as.matrix(tmb_report$cov.fixed[sigma_ind,sigma_ind])
  re <- get_reTrms(formula, data)

  list(names=sigma_names(formula, data), estimates=sigma_est, covMat=sigma_cov,
       invlink="exp", invlinkGrad="exp", n_obs=nrow(data),
       n_levels=lapply(re$flist, function(x) length(levels(x))), cnms=re$cnms,
       levels=rownames(re$Zt))
}

get_fixed_names <- function(tmb_report){
  out <- names(tmb_report$par.fixed)
  if(is.null(out)) out <- 1:length(tmb_report$par.fixed)
  out
}

print_randvar_info <- function(object){
  group_info <- paste0(names(object$n_levels), ", ",
                       unlist(object$n_levels), collapse="; ")

  val <- do.call(object$invlink, list(object$estimates))

  disp <- data.frame(Groups=names(object$cnms), Name=unlist(object$cnms),
                     Variance=round(val^2,3), Std.Dev.=round(val,3))
  cat("Random effects:\n")
  print(disp, row.names=FALSE)
  #below needs to be corrected for missing values at some point
  #cat(paste0("Number of obs: ",object$n_obs,", groups: ",group_info,"\n"))
}

check_no_support <- function(formula_list){
  has_bars <- any(sapply(formula_list, function(x) !is.null(lme4::findbars(x))))
  if(has_bars){
    stop("This function does not support random effects", call.=FALSE)
  }
}

fit_TMB <- function(model, data, params, random,
                    starts, method, ...){

  fixed_sub <- names(params)[!names(params) %in% random]
  nfixed <- length(unlist(params[fixed_sub]))
  list_fixed_only <- params[fixed_sub]
  plengths <- sapply(list_fixed_only, length)
  starts_order <- rep(fixed_sub, plengths)

  if(!is.null(starts)){
    if(length(starts) != nfixed){
      stop(paste("The number of starting values should be", nfixed))
    }
    list_fixed_only <- params[fixed_sub]
    list_fixed_only <- utils::relist(starts, list_fixed_only)
    params <- replace(params, names(list_fixed_only), list_fixed_only)
  }

  tmb_mod <- TMB::MakeADFun(data = c(model = model, data),
                            parameters = params,
                            random = random,
                            silent=TRUE,
                            DLL = "unmarked_TMBExports")
  tmb_mod$starts_order <- starts_order

  opt <- optim(tmb_mod$par, fn=tmb_mod$fn, gr=tmb_mod$gr, method=method, ...)

  sdr <- TMB::sdreport(tmb_mod, getJointPrecision=TRUE)
  sdr$par <- tmb_mod$par

  AIC = 2 * opt$value + 2 * nfixed

  list(opt=opt, TMB=tmb_mod, sdr=sdr, AIC=AIC)
}

get_coef_info <- function(tmb_report, type, names, idx){
  no_sigma <- !grepl("lsigma", get_fixed_names(tmb_report))
  fixed <- tmb_report$par.fixed[no_sigma] #take out sigmas
  fixed <- fixed[idx]
  names(fixed) <- names
  rand <- get_b_vector(tmb_report, type)
  ests <- c(fixed, rand)
  covMat <- get_joint_cov(tmb_report, type)
  list(ests=ests, cov=covMat)
}

setMethod("sigma", "unmarkedEstimate", function(object, level=0.95, ...){
  rinf <- object@randomVarInfo
  if(length(rinf)==0){
    stop("No random effects in this submodel", call.=FALSE)
  }
  z <- qnorm((1-level)/2, lower.tail = FALSE)
  vals <- rinf$estimates
  ses <- sqrt(diag(rinf$covMat))
  lower <- vals - z*ses
  upper <- vals + z*ses
  Groups <- names(rinf$cnms)
  Name <- unlist(rinf$cnms)
  data.frame(Model=object@short.name, Groups=Groups, Name=Name, sigma=exp(vals),
             lower=exp(lower), upper=exp(upper))
})

setMethod("sigma", "unmarkedFit", function(object, type, level=0.95, ...){
  if(!missing(type)){
    return(sigma(object[type], level=level))
  }
  ests <- object@estimates@estimates
  has_rand <- sapply(ests, function(x) length(x@randomVarInfo)>0)
  if(!any(has_rand)){
    stop("No random effects in this model", call.=FALSE)
  }
  ests <- ests[has_rand]

  out_list <- lapply(ests, sigma, level=level)
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
})


setGeneric("randomTerms", function(object, ...) standardGeneric("randomTerms"))


setMethod("randomTerms", "unmarkedEstimate", function(object, level=0.95, ...){

  rv <- object@randomVarInfo
  if(length(rv)==0){
    stop("No random effects in this submodel", call.=FALSE)
  }

  Groups <- lapply(1:length(rv$cnms), function(x){
                  gn <- names(rv$cnms)[x]
                  rep(gn, rv$n_levels[[gn]])
          })
  Groups <- do.call(c, Groups)

  Name <-  lapply(1:length(rv$cnms), function(x){
                  gn <- names(rv$cnms)[x]
                  var <- rv$cnms[[x]]
                  rep(var, rv$n_levels[[gn]])
          })
  Name <- do.call(c, Name)

  rv_idx <- !1:length(object@estimates) %in% object@fixed
  b_var <- object@estimates[rv_idx]
  b_se <- sqrt(diag(object@covMat[rv_idx,rv_idx,drop=FALSE]))

  z <- qnorm((1-level)/2, lower.tail = FALSE)
  lower <- b_var - z*b_se
  upper <- b_var + z*b_se

  data.frame(Model=object@short.name, Groups=Groups, Name=Name, Level=rv$levels,
                    Estimate=b_var, SE=b_se, lower=lower, upper=upper)
})


setMethod("randomTerms", "unmarkedFit", function(object, type, level=0.95, ...){

  if(!missing(type)){
    return(randomTerms(object[type], level))
  }

  has_random <- sapply(object@estimates@estimates,
                       function(x) length(x@randomVarInfo) > 0)
  if(!any(has_random)){
    stop("No random effects in this model", call.=FALSE)
  }
  keep <- object@estimates@estimates[has_random]
  out <- lapply(keep, randomTerms, level=level)
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
})

get_ranef_inputs <- function(forms, datalist, dms, Zs){
  stopifnot(!is.null(names(datalist)))
  mods <- names(datalist)
  ngv <- lapply(forms, get_group_vars)
  names(ngv) <- paste0("n_group_vars_",mods)
  ngroup <- mapply(get_nrandom, forms, datalist, SIMPLIFY=FALSE)
  names(ngroup) <- paste0("n_grouplevels_",mods)
  names(dms) <- paste0("X_", mods)
  names(Zs) <- paste0("Z_", mods)

  dat <- c(ngv, ngroup, dms, Zs)

  beta <- lapply(dms, function(x) rep(0, ncol(x)))
  names(beta) <- paste0("beta_", mods)
  b <- lapply(ngroup, function(x) rep(0, sum(x)))
  names(b) <- paste0("b_", mods)
  lsigma <- lapply(ngv, function(x) rep(0, x))
  names(lsigma) <- paste0("lsigma_", mods)

  pars <- c(beta, b, lsigma)

  rand_ef <- paste0(names(b))[sapply(forms, has_random)]
  if(length(rand_ef) == 0) rand_ef <- NULL

  list(data=dat, pars=pars, rand_ef=rand_ef)
}

add_covariates <- function(covs_long, covs_short, n){

  if(is.null(covs_short)){
    return(covs_long)
  }

  if(is.null(covs_long)){
    covs_long <- data.frame(.dummy=rep(1, n))
  } else {
    stopifnot(nrow(covs_long) == n)
  }

  exp_factor <- nrow(covs_long) / nrow(covs_short)
  stopifnot(exp_factor > 1)

  rep_idx <- rep(1:nrow(covs_short), each=exp_factor)

  to_add <- covs_short[rep_idx, ]
  stopifnot(nrow(covs_long) == nrow(to_add))

  cbind(covs_long, to_add)
}

vcov_TMB <- function(object, type, fixedOnly){

  if(!missing(type)){
    return(vcov(object[type], fixedOnly=fixedOnly))
  }

  v <- get_joint_cov(TMB::sdreport(object@TMB, getJointPrecision=TRUE))
  no_sig <- !grepl("lsigma_",colnames(v))
  v <- v[no_sig, no_sig]
  colnames(v) <- rownames(v) <- names(coef(object, fixedOnly=FALSE))

  if(fixedOnly){
    no_re <- !grepl("b_", colnames(v))
    v <- v[no_re, no_re]
  }
  v
}
