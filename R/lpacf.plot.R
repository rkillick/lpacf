lpacf.plot=function(lpacf,atTime=NULL,atLag=NULL,SaveToFile=FALSE,alpha=0.95,...){
	# plots the lpacf at row atTime and/or column atLag
	# both atTime and atLag can be vectors

	# check inputs
  if(class(lpacf)!='lpacf'){stop('lpacf must be of class lpacf')}
  if((alpha<=0)|(alpha>=1)){stop('alpha must be strictly between 0 and 1.')}
  ci=qnorm((1 + alpha)/2)/sqrt(lpacf$binwidth) # only valid if allpoints=F
  if(any(names(lpacf)=='the.x1')){warning("Confidence intervals for the ends of the data are incorrect.  RK needs to code this.")}
  
	if(ncol(lpacf$lpacf)>nrow(lpacf$lpacf)){stop('cannot have more lags than time in lpacf, must be a matrix of dimension n x lag')}
	if(any(atTime<=0)){stop('atTime cannot contain 0 or negative elements')}
	if(any(atTime>nrow(lpacf$lpacf))){stop('atTime cannot contain elements greater than the number of rows in lpacf')}
	if(any(atLag<=0)){stop('atLag cannot contain 0 or negative elements')}
	if(any(atLag>ncol(lpacf$lpacf))){stop('atLag cannot contain elements greater than the number of columns in lpacf')}

	# now do the plots
	if(SaveToFile){pdf(file='Rplotlpacf%03d.pdf')}

	for(i in atTime){
		ts.plot(lpacf$lpacf[i,],main=paste('lpacf at time',lpacf$the.x[i]),xlab='Lag',ylab='Local Partial ACF',ylim=c(-1,1),...)
    abline(h=c(-ci,ci),lty=2,col='blue')
		if(!SaveToFile){ readline("Pause. Press <Enter> to continue...")}
	}
	for(i in atLag){
		ts.plot(lpacf$lpacf[,i],main=paste('lpacf at lag',i),xlab='Time',ylab='Local Partial ACF',ylim=c(-1,1),...)
		abline(h=c(-ci,ci),lty=2,col='blue')
		if(!SaveToFile){ readline("Pause. Press <Enter> to continue...")}
	}
	if(SaveToFile){dev.off()}

}
