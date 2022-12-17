# GEV functions
library(evd)


###
# log posterior of GEV function
###
gev_log_post = function(theta,
									 Y, W,
									 mu, inv_Sigma){
	
	# likelihood
	log_post = sum(dgev(Y, W %*% theta[1:2], exp(theta[3]),theta[4],log=TRUE))

	# prior
	log_post = log_post - 0.5*t(theta-mu) %*% inv_Sigma %*% (theta-mu)
	
	return(log_post)
}
