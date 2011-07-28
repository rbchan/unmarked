sim.data<-function(N=30,T=5,sigma=.6,lam0=1.5,knownID=0,lc=10,buff=3,sslower=0,ssupper=15){

  #  lc <- 10
    coords <- seq(sslower+buff, ssupper-buff, length=lc)
    X <- cbind(x=rep(coords, each=lc), y=rep(coords, times=lc))
 ng<- npts<-nrow(X)
    # Home range centers
    sx <- runif(N, sslower,ssupper)
    sy <- runif(N, sslower,ssupper)
    S <- cbind(sx, sy)
    D <- e2dist1(S, X)
# plot(X,pch=20)
# points(S,pch=20,col="red")
 lam<- lam0*exp(-(D*D)/(2*sigma*sigma))
 zlst<-list()
 y<-matrix(NA,nrow=npts,ncol=T)
 for(i in 1:T){
 z<-rpois(prod(dim(D)),lam[1:length(lam)])
 z<-matrix(z,nrow=nrow(lam),ncol=ncol(lam),byrow=FALSE)
 zlst[[i]]<-z
 }
ncaps<-matrix(NA,nrow=N,ncol=T)
for(i in 1:T){
ncaps[,i]<-apply(zlst[[i]],1,sum)
}
capind<- (1:N)[apply(ncaps,1,sum)>0]

if(knownID>0){
if(length(capind)<=knownID){
 cat("ERROR of a serious kind",fill=TRUE)
 gooseegg<-NULL
 return(gooseegg)
}
known.kp<-capind[1:knownID]
unknown.kp<- is.na(match(1:N,known.kp))
zknown<-list()
for(i in 1:T){
 zknown[[i]]<-zlst[[i]][known.kp,]  # knownID x traps matrix
 y[,i]<-t(zlst[[i]][unknown.kp,])%*%rep(1,N-knownID)
 }
}
if(knownID==0){
zknown=NULL
for(i in 1:T){
 y[,i]<-t(zlst[[i]])%*%rep(1,N)
}
}

if(knownID>0){
zarr<-array(0,c(knownID,ng,T))
 for(i in 1:T){
 zarr[,,i]<- zknown[[i]]
  }
}
else{
zarr<-NULL
}



list(y=y,N=N,T=T,sigma=sigma,lam0=lam0,X=X,S=S,knownID=knownID,zknown=zknown,zarr=zarr)


}
