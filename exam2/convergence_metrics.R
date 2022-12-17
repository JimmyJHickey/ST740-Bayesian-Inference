library(coda)


calc_ESS = function(param_trace, max_iteration, burnin){
	
	n = dim(param_trace)[1]
	m = dim(param_trace)[2]
	
	
	ESS = array(0,c(n,m))
	
	for(i in 1:n){
		samps = mcmc(t(param_trace[i,,burnin:max_iteration]), thin=1)
		ESS[i,] = effectiveSize(samps)
	}
	
	return(ESS)
}

calc_geweke= function(param_trace, max_iteration, burnin){
	
	n = dim(param_trace)[1]
	m = dim(param_trace)[2]
	
	library(coda)
	geweke = array(0,c(n,m))
	
	for(i in 1:n){
		samps <- mcmc(t(param_trace[i,,burnin:max_iteration]), thin=1)
		geweke[i,] = geweke.diag(samps)$z
	}
	
	return(geweke)
}