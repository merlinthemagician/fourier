fourierNames=list.files(path="~/c/fourier/results",pattern="*Fourier.dat");

fourierData=lapply(fourierNames,function(x) read.table(paste("~/c/fourier/results/",x, sep=""),sep="\t", colClasses="complex"));
