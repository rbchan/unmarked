setGeneric("posteriorSamples", function(object, nsims, ...){
             standardGeneric("posteriorSamples")
          })

setClass("unmarkedPostSamples",
         representation(numSites="numeric",
                        numPrimary="numeric",
                        nsims="numeric",
                        samples="array")
         )

setMethod("posteriorSamples", "unmarkedRanef", function(object, nsims=100, ...)
{

  N <- dim(object@post)[1]
  K <- dim(object@post)[2]
  T <- dim(object@post)[3]

  out <- array(NA, c(N, T, nsims))

  for (n in 1:N){
    for (t in 1:T){
        out[n, t, ] <- sample(0:(K-1), nsims, replace=TRUE,
                              prob=object@post[n,,t])
    }
  }
  new("unmarkedPostSamples", numSites=N, numPrimary=T, nsims=nsims,
      samples=out)

})

setMethod("posteriorSamples", "unmarkedFit", function(object, nsims=100, ...)
{
  ran <- ranef(object)
  posteriorSamples(ran, nsims)
})

setMethod("show", "unmarkedPostSamples", function(object)
{

  #tdim <- character(0)
  #if(object@numPrimary>1){
  tdim <- paste0("x ", object@numPrimary, " primary periods")
  #}

  cat("Posterior samples from unmarked model\n")
  cat(paste(object@numSites, "sites", tdim, "x", object@nsims, "sims\n"))
  cat(paste0("Showing first 5 sites and first 3 simulations\n",
      "To see all samples, use print()\n"))
  
  print(object@samples[1:5,,1:3])

})

print.unmarkedPostSamples <- function(x, ...){
  print(x@samples)
}

setMethod("[", c("unmarkedPostSamples","ANY","ANY","ANY"), 
          function(x, i, j, k, drop = FALSE)
{
  x@samples[i,j,k, drop = drop]
})
