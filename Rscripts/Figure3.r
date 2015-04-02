## Figure 3:
loadNoise=function(pth, pat) {
    names=list.files(path=pth,pattern=pat);
    return(lapply(names,function(x) read.table(paste(pth,x, sep="/"),sep=" ", colClasses="numeric")))
}

NoiseList=loadNoise("~/c/fourier/results/w5/", "*t0000.dat");

listAllij=function(Nlist, i, j) {
    return(sapply(Nlist, function (x) x[i+1,j+1]))
}

R=0:20

plot(R, sapply(R, function(x) cov(listAllij(NoiseList,0,0),listAllij(NoiseList,0,x))/var(listAllij(NoiseList,0,0))))
