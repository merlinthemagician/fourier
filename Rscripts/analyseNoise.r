fourierNamesW5=list.files(path="~/c/fourier/results/w5",pattern="*Fourier.dat");

fourierDataW5=lapply(fourierNames,function(x) read.table(paste("~/c/fourier/results/w5/",x, sep=""),sep="\t", colClasses="complex"));

noiseNames=list.files(path="~/c/fourier/results/w5",pattern="*[0-9].dat");
noiseData=lapply(noiseNames,function(x) read.table(paste("~/c/fourier/results/w5/",x, sep=""),sep=" ", colClasses="numeric"));

Sstat=function(eps, tau,lambda, k) {
    return(eps/(tau*(1+lambda*lambda*k*k)))
}

cMuNu=function(lambda, hx, N, mu, nu) {
    return(1-2*lambda^2/(hx^2) *(cospi(2*mu/N) + cospi(2*nu/N))-2)
}

SstatDisc=function(eps,tau,lambda, hx, N, mu, nu) {
    return(eps/(tau*cMuNu(lambda,hx,N,mu,nu)))    
}

4096/var(sapply(fourierDataW5[c(1:2500)],function(x) x[1,1]))

