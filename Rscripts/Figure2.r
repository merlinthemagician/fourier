## Figure 2:
loadComplex=function(pth, pat) {
    names=list.files(path=pth,pattern=pat);
    return(lapply(names,function(x) read.table(paste(pth,x, sep="/"),sep="\t", colClasses="complex")))
}

cMuNu=function(lambda, hx, N, mu, nu) {
    return(1-2*lambda^2/(hx^2) *(cospi(2*mu/N) + cospi(2*nu/N)-2))
}

timeCorrScl=function(eps,tau, lambda, hx, N, mu, nu, diffT) {
    cMN=cMuNu(lambda, hx, N, mu, nu)
    return(eps/(tau*cMN)*exp(-cMN/tau*diffT))
}

Klist=loadComplex("~/c/fourier/results/", "*t*.dat")

fudgeComplexCovariance=function(K1,K2,sel) {
    return(Mod(Reduce('+', (K1[[sel]]-mean(K1[[sel]]))*Conj(K2[[sel]]-mean(K2[[sel]]))))/(length(K1[[sel]]-1)))
}


times=0:(length(Klist)-1)/2

eps=5
tau=1
lambda=3
hx=1
N=64

## k00
plot(times, sapply(Klist, function(x) cov(Re(Klist[[1]]$V1),Re(x$V1))/4096))
plot(function(x) timeCorrScl(eps,tau, lambda, hx, N, 0, 0, x),0,times[length(times)], add=T)

## k02
dev.new()
plot(times, sapply(Klist, function(x) fudgeComplexCovariance(Klist[[1]],x,"V2")/4096))

plot(function(x) timeCorrScl(eps,tau, lambda, hx, N, 0, 2, x),0,times[length(times)], add=T)

## k04
dev.new()
plot(times, sapply(Klist, function(x) fudgeComplexCovariance(Klist[[1]],x,"V3")/4096))

plot(function(x) timeCorrScl(eps,tau, lambda, hx, N, 0, 4, x),0,times[length(times)], add=T)


## k06
dev.new()
plot(times, sapply(Klist, function(x) fudgeComplexCovariance(Klist[[1]],x,"V4")/4096))

plot(function(x) timeCorrScl(eps,tau, lambda, hx, N, 0, 6, x),0,times[length(times)], add=T)


## Combined
## dev.new()

## plot(times, list(sapply(Klist, function(x) cov(Re(Klist[[1]]$V1),Re(x$V1))/4096),sapply(Klist, function(x) fudgeComplexCovariance(Klist[[1]],x,"V2")/4096)))

## plot(function(x) timeCorrScl(eps,tau, lambda, hx, N, 0, 0, x),0,times[length(times)], add=T)

## plot(function(x) timeCorrScl(eps,tau, lambda, hx, N, 0, 2, x),0,times[length(times)], add=T)
