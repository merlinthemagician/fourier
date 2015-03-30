fourierNames=list.files(path="~/c/fourier/results",pattern="*Fourier.dat");

fourierData=lapply(fourierNames,function(x) read.table(paste("~/c/fourier/results/",x, sep=""),sep="\t", colClasses="complex"));

noiseNames=list.files(path="~/c/fourier/results",pattern="*[0-9].dat");
noiseData=lapply(noiseNames,function(x) read.table(paste("~/c/fourier/results/",x, sep=""),sep="\t"));
