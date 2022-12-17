source("update_params.R")
source("gev_log_post.R")
# GEV functions
library(evd)


MCMC = function( # data
	Y,
	W,
	# mcmc parameters
	max_iteration = 10000,
	burnin = 5000,
	# initial parameters with defaults
	mu0 = rep(0,4),
	Sigma0 = 1000*diag(4),
	# hyper parameters with defaults
	c = 100,
	Psi = diag(4),
	nu = 3
){
	start_time = proc.time()
	n = dim(Y)[1]
	m = dim(Y)[2]
	
	###
	# initial parameters
	###
	mu = mu0
	Sigma = Sigma0
	
	# initialize based on MLE 
	theta = matrix(0, n, 4)
	p = dim(theta)[2]
	current_log_post = rep(0, n)
	for( loc in 1:n ){
		theta_mle = fgev(Y[loc,],nsloc=W[,2],warn.inf=FALSE,std.err = FALSE)
		# theta_mle = fgev(Y[loc,],nsloc=W[,2],warn.inf=FALSE)
		theta[loc,] = theta_mle$estimate
		theta[loc,3] = log(theta[loc,3]) # transform sigma to log scale
		current_log_post[loc] = gev_log_post(theta[loc,],Y[loc,], W, mu, solve(Sigma))
	}
	
	# start trace matrices
	mu_trace = array(0, c(1, 4, max_iteration))
	mu_trace[,,1] = mu
	
	Sigma_trace = array(0, c(4,4, max_iteration))
	Sigma_trace[,,1] = Sigma
	
	
	theta_trace=array(0,c(n,4,max_iteration))
	theta_trace[,,1] = theta
	
	
	# MH tuning parameters
	# acceptance ratio lower and upper bound
	ar_lower = 0.3 
	ar_upper = 0.5
	# acceptance ratio lower and upper multipliers
	ar_lower_mult = 0.8
	ar_upper_mult = 1.2
	
	# number of iterations before tuning MH ratio
	tuning_iterations = 50
	
	attempted = accepted = 0
	theta_mle = fgev(Y[1,],nsloc=W[,2],warn.inf=FALSE)
	MH = theta_mle$std.err            # Standard errors as initial guesses
	MH[3] = MH[3]/theta_mle$estimate[3]  # Delta method for log(scale)
	MH = diag(MH)
	
	###
	# run MCMC
	###
	for(iteration in 2:max_iteration){
		
		if(iteration %% 100 == 0){print(iteration)}
		
		inv_Sigma = solve(Sigma)
		
		###
		# Metropolis-Hastings to update theta
		###
		for(loc in 1:n){
			
			
			attempted = attempted + 1
			
			candidate = theta[loc,]
			candidate = mvrnorm(1,theta[loc, ], MH)
			candidate_log_post = gev_log_post(candidate, Y[loc,], W, mu, inv_Sigma)
			
			# accept or reject new candidate theta vector
			if( !is.na(candidate_log_post) &
					log(runif(1)) < candidate_log_post - current_log_post[loc] ){
				
				accepted = accepted + 1
				theta[loc, ] = candidate
				current_log_post[loc] = candidate_log_post
			}
		}
		
		
		
		
		# tune MH step size
		if( iteration < burnin){
			
			if(attempted > tuning_iterations){
				ar = accepted / attempted
				if(ar < ar_lower){MH = MH*ar_lower_mult}
				if(ar > ar_upper){MH = MH*ar_upper_mult}
				
				accepted = attempted = 0
			}
		}
		
		
		###
		# Gibbs sampling for mu, Sigma
		###
		mu = update_mu(inv_Sigma, theta, c)
		Sigma = update_Sigma(theta, mu, Psi, nu)
		
		# save current value to trace
		mu_trace[,,iteration] = mu
		Sigma_trace[,,iteration] = Sigma
		theta_trace[,,iteration] = theta
	}
	
	end_time = proc.time()
	return(list(mu_trace = mu_trace,
							Sigma_trace = Sigma_trace,
							theta_trace = theta_trace,
							acceptance_rate = accepted / attempted,
							time=end_time - start_time))
	
}

