# sample from inverse gamma
library(MCMCpack)
# sample from dirichlet
library(LaplacesDemon)

# draw from full conditional of theta_i
update_theta_i = function(Yi,
													Mi, alpha){
	
	return(MCMCpack::rdirichlet(1, exp(Mi) * alpha + Yi))
}


# draw from full conditional of mu
update_mu = function(M, sigmasq,
										 m0, ssq){
	
	n = length(M)
	variance = 1/(n / sigmasq + 1 / ssq)

	return(rnorm(1, variance * (sum(M)/sigmasq + m0/ssq), variance))
	
}


# draw from full conditional of sigmasq
update_sigmasq = function(M, mu,
													a, b){
	
	n = length(M)
	sumsq = sum((M - mu)^2)
	
	return(MCMCpack::rinvgamma(1, n/2 + a, sumsq/2 + b))
}


# Mi log post
m_i_log_post= function(m_i, theta_i, alpha,
												mu, sigmasq){
	log_post = dnorm(m_i, mu, sqrt(sigmasq), log = TRUE) + 
							LaplacesDemon::ddirichlet(theta_i, exp(m_i) * alpha, log = TRUE)
	
	return(log_post)
}
