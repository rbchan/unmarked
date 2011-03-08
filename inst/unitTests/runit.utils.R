test.formatLong <- function() {
  df <- read.csv(system.file("csv","frog2001pcru.csv", package = "unmarked"))
  umf <- formatLong(df, type = "unmarkedFrameOccu")
  ## Add some assertions...
}





test.rowProds <- function()
{
    m <- matrix(1:9/10, 3)
    rowProds <- unmarked:::rowProds
    checkEquals(rowProds(m), apply(m, 1, prod))    
}





test.tranProbs <- function()
{
    N <- 0:100
    omega <- 0.5
    gamma <- 2
    delta <- 5
    tranProbs <- unmarked:::tranProbs
    tranProbsR <- unmarked:::tranProbsR
    tp1 <- tranProbs(N, omega, gamma, delta, "constant")
    tp1r <- tranProbsR(N, omega, gamma, delta, "constant")
    checkEquals(tp1, tp1r)
    
    checkTrue(all.equal(colSums(tp1), rep(1, 101)))


    tp1 <- tranProbs(Nr=0:4, omegaR=0.5, gammaR=1, deltaR=1, 
        dynamicsR="constant")
    tp1r <- tranProbsR(0:4, omega=0.5, gamma=1, delta=1, 
        dynamics="constant")
    checkEquals(tp1, tp1r)
    
    tp2 <- tranProbs(Nr=0:4, omegaR=0.5, gammaR=1, deltaR=3, 
        dynamicsR="constant")
    tp2r <- tranProbsR(0:4, omega=0.5, gamma=1, delta=3, 
        dynamics="constant")
    checkEquals(tp2, tp2r)
    
    tp3 <- tranProbs(Nr=0:4, omegaR=0.5, gammaR=1, deltaR=1, 
        dynamicsR="autoreg")
    tp3r <- tranProbsR(0:4, omega=0.5, gamma=1, delta=1, 
        dynamics="autoreg")
    checkEquals(tp3, tp3r)
    
    tp4 <- tranProbs(Nr=0:4, omegaR=0.5, gammaR=1, deltaR=5, 
        dynamicsR="autoreg")
    tp4r <- tranProbsR(0:4, omega=0.5, gamma=1, delta=5, 
        dynamics="autoreg")
    checkEquals(tp4, tp4r)
    

    
}
    
    
