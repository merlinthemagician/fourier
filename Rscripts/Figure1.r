fourierNamesW5t0=list.files(path="~/c/fourier/results/w5",pattern="*t0000Fourier.dat");

fourierDataW5t0=lapply(fourierNamesW5t0,function(x) read.table(paste("~/c/fourier/results/w5/",x, sep=""),sep="\t", colClasses="complex"));

## noiseNames=list.files(path="~/c/fourier/results/w5",pattern="*[0-9].dat");
## noiseData=lapply(noiseNames,function(x) read.table(paste("~/c/fourier/results/w5/",x, sep=""),sep=" ", colClasses="numeric"));

Sstat=function(eps, tau,lambda, k) {
    return(eps/(tau*(1+lambda*lambda*k*k)))
}

cMuNu=function(lambda, hx, N, mu, nu) {
    return(1-2*lambda^2/(hx^2) *(cospi(2*mu/N) + cospi(2*nu/N)-2))
}

cMuNuCont=function(lambda, k) {
    return(1+lambda*lambda*k*k)
}

SstatDisc=function(eps,tau,lambda, hx, N, mu, nu) {
    return(eps/(tau*cMuNu(lambda,hx,N,mu,nu)))    
}

listAllMuNu=function(Klist, mu, nu) {
    return(sapply(Klist, function (x) x[mu+1,nu+1]))
}

compCov=function(Klist, mu, nu) {
    AllKMuNu=listAllMuNu(Klist, mu, nu)
    EKMuNu=mean(AllKMuNu)
    return(Mod(Reduce('+', (AllKMuNu-EKMuNu)*Conj(AllKMuNu-EKMuNu)))/(length(Klist)-1))
}

## Plot theoretical results
## Discrete:
##plot(8*pi*pi/4096*K20*K20,sapply(K20,function(x) 1/SstatDisc(5,1,3,1,64,x,x)),ylim=c(0,15))
## Continuous:
##plot(function(x) cMuNuCont(3,sqrt(x))/5,0,8, add=T)

##dev.new()
plot(8*pi*pi/4096*K20*K20,sapply(K20,function(x) compCov(fourierDataW5t0, x, x) ),ylim=c(0,15))

plot(8*pi*pi/4096*K20*K20,sapply(K20,function(x) 4096/compCov(fourierDataW5t0, x, x) ),ylim=c(0,15),col="red")

par(new=TRUE)

plot(8*pi*pi/4096*K20*K20,sapply(K20,function(x) 1/SstatDisc(5,1,3,1,64,x,x)),ylim=c(0,15), axes=FALSE)

plot(function(x) cMuNuCont(3,sqrt(x))/5,0,8, add=T)

## loadFourier=function(pth, pat) {
##     names=list.files(path=pth,pattern=pat);

##     return(lapply(names,function(x) read.table(paste(pth,x, sep="/"),sep="\t", colClasses="complex")))
## }

