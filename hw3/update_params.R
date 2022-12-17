# sample from multivariate normal
library(MASS)

# sample from inverse gamma
library(MCMCpack)

#sample from inverse wishart
library(LaplacesDemon)

# draw from the full conditional of mu
update_mu = function(Sigma, beta, c){
	
	inv_Sigma = solve(Sigma)
	n = dim(beta)[1]
	
	cov_mat = solve(n * inv_Sigma + 1/c * diag(2))
	
	mmm = inv_Sigma %*% colSums(beta)
	
	return( mvrnorm(1, cov_mat %*% mmm, cov_mat))
}


# draw from full conditional of sigmasq
update_sigmasq = function(Z, W, beta, a, b){
	n = dim(Z)[1]
	m = dim(Z)[2]
	
	iga = a + n*m/2
	
	igb = 0
	for(s in 1:n){
		temp = Z[s,] - W %*% beta[s,]
		igb = igb +(  t(temp) %*% temp  )
	}
	
	return(rinvgamma(1, iga, b + 1/2 * igb))
	
}


# draw from full conditional of Sigma
update_Sigma = function(beta, mu, Psi, nu){
	
	
	n = dim(beta)[1]
	
	S = matrix(0, 2,2)
	
	for(s in 1:n){
		temp = beta[s, ]-mu
		S = S + temp %*% t(temp)
	}
	
	
	
	return(rinvwishart(nu + n, Psi + S))
}


update_beta = function(Zs, W, mu, Sigma, sigmasq){
	inv_Sigma = solve(Sigma)
	cov_mat = solve(1/sigmasq * t(W) %*% W + inv_Sigma )
	mean_mat = 1/sigmasq * t(W) %*% Zs + inv_Sigma %*% mu
	
	return(mvrnorm(1, cov_mat %*% mean_mat, cov_mat))
}


# 
# m = 4
# n = 3
# 
# Z_test = rbind(c(1,1,1), c(2,2,2), c(3,3,3), c(4,4,4))
# 
# X_test = c(1,1,1)
# W_test = matrix(c(rep(1, n), X_test), n, 2)
# a=0.5
# b=0.5
# 
# Sigma = diag(2)
# 
# beta = rbind(c(10,10), c(20,20), c(30, 30), c(40, 40))
# 
# 
# update_simgasq(Z_test, W_test, beta, a, b)
# 
# Psi = diag(2)
# nu = 3
# 
# mu = c(1,2)
# 
# update_Sigma(beta, mu, Psi, nu)
# 
# sigmasq = 2

# solve(1/sigmasq * t(W_test) %*% W_test + solve(Sigma) ) %*% t(W_test) %*% Z_test[1,]


# update_beta(Z_test[1,], W_test, mu, Sigma, sigmasq)

# Sigma = diag(2)
# 
# realn = 3
# beta = t(matrix(c(10, 3,
# 							6, 15,
# 							44, 1), 2, realn))
# c = 0.5
# 
# 
# Sigmainv = solve(Sigma)
# 
# update_mu(solve(Sigma), beta, c)
