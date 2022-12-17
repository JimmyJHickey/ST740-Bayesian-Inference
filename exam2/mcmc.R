
source("update_params.R")
library(MCMCpack)
library(stats)

log_like <- function(Y, theta){
	
	n = dim(Y)[2]
	ll = rep(0, n)
	
	
	for(i in 1:n){
		ni = sum(Y[,i])
		ll[i] = dmultinom(Y[,i], ni, theta[i,], log=TRUE)
	} 
	
	return(ll)
}


MCMC = function(Y,
								max_iteration = 2000,
								burnin = 200,
								# if M_fixed is TRUE then M will be set to M0
								# otherwise its distribution will be estimated
								est_M = TRUE,
								M0 = 0,
								MH0 = 0.1,
								# hyper parameters with defaults
								m0 = 50,
								ssq = 1,
								a = 1,
								b = 1,
								seed = 1978){
	
	start_time = proc.time()
	
	set.seed(seed)


	
	n = dim(Y)[2]
	m = dim(Y)[1]
	
	# sample proportions of each taxa
	# alpha = (rowSums(Y+1)) / (sum(colSums(Y+1)))
	alpha = rowSums(Y) / sum(Y)
	
	
	mu_trace = rep(0, max_iteration)
	mu_trace[1] = rnorm(1, m0, ssq)
	
	
	sigmasq_trace = rep(0, max_iteration)
	sigmasq_trace[1] = rinvgamma(1,a,b)
	
	
	M_trace = array(0, c(n,max_iteration))
	
	
	MH = rep(MH0, n)
	accepted = attempted = rep(0,n)
	
	# MH tuning parameters
	# acceptance ratio lower and upper bound
	ar_lower = 0.3 
	ar_upper = 0.5
	# acceptance ratio lower and upper multipliers
	ar_lower_mult = 0.8
	ar_upper_mult = 1.2
	
	# number of iterations before tuning MH ratio
	tuning_iterations = 20
	
	
	if(est_M){
		M_trace[,1] = rnorm(n, mu_trace[1], sigmasq_trace[1])
	} else {
		M_trace[,1] = M0
	}
	
	theta_trace = array(0,c(n,m,max_iteration))
	for(i in 1:n){
		theta_trace[i,,1] = rdirichlet(1, exp(M_trace[i,1]) * alpha)
	}
	
	
	logpdf1 = logpdf2 = 0
	
	# mcmc
	for(iteration in 2:max_iteration){
		
		if(iteration %% 100 == 0) {print(paste("iteration ", iteration))}
		
		# Gibbs for thetai, mu, sigmasq
		# Metropolis for Mi


		if(est_M){
			mu_trace[iteration] = update_mu(M_trace[,iteration-1], sigmasq_trace[iteration-1], m0, ssq)
			
			sigmasq_trace[iteration] = update_sigmasq(M_trace[,iteration-1], mu_trace[iteration], a, b)
			
			# update Mi
			for(i in 1:n){
			
				attempted[i] = attempted[i] + 1
				
				current_Mi = M_trace[i,iteration-1]
				M_trace[i,iteration] = current_Mi

				current_log_post = m_i_log_post(current_Mi,
																				theta_trace[i,,iteration-1],
																				alpha,
																				mu_trace[iteration],
																				sigmasq_trace[iteration])
				
				candidate_Mi = rnorm(1, current_Mi, MH[i])
				candidate_log_post = m_i_log_post(candidate_Mi,
																					theta_trace[i,,iteration-1],
																					alpha,
																					mu_trace[iteration],
																					sigmasq_trace[iteration])

				
				# accept or reject new candidate theta vector
				if( !is.na(candidate_log_post) &
						log(runif(1)) < candidate_log_post - current_log_post ){
				
					accepted[i] = accepted[i] + 1
					M_trace[i, iteration] = candidate_Mi
				}
				
				
				# tune MH step size
				if( iteration < burnin & attempted[i] > tuning_iterations){
					ar = accepted[i] / attempted[i]
					if(ar < ar_lower){MH[i] = MH[i]*ar_lower_mult}
					if(ar > ar_upper){MH[i] = MH[i]*ar_upper_mult}
					
					accepted[i] = attempted[i] = 0
				}
					
			} # for 1:n
		} # if M
		
		
		
		for(i in 1:n){
			theta_trace[i,,iteration] = update_theta_i(Y[,i], M_trace[i,iteration], alpha)
		}
		
		
		if(iteration > burnin){
			# WAIC calculations
			ll = log_like(Y, theta_trace[,,iteration])
			logpdf1 = logpdf1 + (ll)/(max_iteration-burnin)
			logpdf2 = logpdf2 + (ll^2)/(max_iteration-burnin)
		}
		
	} # for MCMC iteration
		
	
	mn_logpdf  <- logpdf1
	var_logpdf <- logpdf2 - logpdf1^2
	pW         <- sum(var_logpdf)
	WAIC       <- list(WAIC=-2*sum(mn_logpdf)+2*pW,pW=pW)
	
	
	
	end_time = proc.time()
	
	
	
	
	return(list(mu_trace = mu_trace,
							sigmasq_trace = sigmasq_trace,
							theta_trace = theta_trace,
							M_trace = M_trace,
							WAIC = WAIC,
							time = end_time - start_time))
		
}
	