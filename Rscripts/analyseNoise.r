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

refind=function(i,n) {
    if(i==1) {return(i);}
    else {return(n-i);}
}

covMatrix=function(M) {
    A=matrix(c(1:n*n),n,n);
    n=dim(M)[1]
    for(i in 1:n) {
        for(j in 1:i) {
            A[i,j]=(sqrt(M[i,j]*M[refind(i,n+2),refind(j,n+2)]));
            A[j,i]=A[i,j]
        }
    }
    return(A);
}

allcov=lapply(fourierDataW5[c(1:2500)],function(x) covMatrix(x))
Reduce('+',sapply(fourierDataW5[c(1:2500)], function(x) (x[10,10]*x[56,56])))/2500

dat2CovMatrix=function(M) {
    return(sapply(M,function(x) Mod(x)*Mod(x)))
}

listdat2CovMatrix=function(listM) {
    return(lapply(listM, dat2CovMatrix))
}

listDat2Cov=function(listM) {
    return(Reduce('+', listdat2CovMatrix(listM))/length(listM))
}

loadFourier=function(pth, pat) {
    names=list.files(path=pth,pattern=pat);

    return(lapply(names,function(x) read.table(paste(pth,x, sep="/"),sep="\t", colClasses="complex")))
}

