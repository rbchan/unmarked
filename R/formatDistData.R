formatDistData <- function (distData, distCol, transectNameCol, dist.breaks, occasionCol,effortMatrix) 
{
  if (!is.numeric(distData[, distCol])) 
    stop("The distances must be numeric")
  transects <- distData[, transectNameCol]
  if (!is.factor(transects)) {
    transects <- as.factor(transects)
    warning("The transects were converted to a factor")
  }
  if (missing(occasionCol)) {
    T <- 1
    occasions <- factor(rep(1, nrow(distData)))
  }
  if (missing(effortMatrix)) {
    effortMatrix <- matrix(nrow=length(unique(transects)),ncol=T,1)
  }
  else {
    occasions <- distData[, occasionCol]
    if (!is.factor(occasions)) {
      occasions <- as.factor(occasions)
      warning("The occasions were converted to a factor")
    }
    T <- nlevels(occasions)
  }
  M <- nlevels(transects)
  J <- length(dist.breaks) - 1
  dist.classes <- levels(cut(distData[, distCol], dist.breaks, 
                             include.lowest = TRUE))
  ya <- array(NA, c(M, J, T), dimnames = list(levels(transects), 
                                              dist.classes, paste("rep", 1:T, sep = "")))
  transect.levels <- levels(transects)
  occasion.levels <- levels(occasions)
  for (i in 1:M) {
    for (t in 1:T) {
      sub <- distData[transects == transect.levels[i] & 
                        occasions == occasion.levels[t], , drop = FALSE]
      ya[i, , t] <- table(cut(sub[, distCol], dist.breaks, 
                              include.lowest = TRUE))
    }
  }
  y <- matrix(ya, nrow = M, ncol = J * T)
  # takes into account the effortMatrix to allow for the insertion of NAs instead of 0s for surveys which were not completed  
  ee <- array(NA, c(M,length(occasion.levels)*(length(dist.breaks)-1)))
  for(i in 1:length(occasion.levels)){
  ee[,((ncol(ee)/length(occasion.levels)*(i-1)+1):(ncol(ee)/length(occasion.levels)*i))] <- matrix(effortMatrix[,i], ncol=(length(dist.breaks)-1), nrow=nrow(effortMatrix))
  }
  ee[ee==0] <- NA
  y <- y * ee
  dn <- dimnames(ya)
  rownames(y) <- dn[[1]]
  if (T == 1) 
    colnames(y) <- dn[[2]]
  else colnames(y) <- paste(rep(dn[[2]], times = T), rep(1:T, each = J), sep = "")
  return(y)
}
