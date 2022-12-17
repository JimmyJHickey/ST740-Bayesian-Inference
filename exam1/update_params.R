# sample from multivariate normal
library(MASS)

#sample from inverse wishart
library(LaplacesDemon)

# draw from the full conditional of mu
update_mu = function(inv_Sigma, theta, c){
	
	n = dim(theta)[1]
	
	cov_mat = solve(n * inv_Sigma + 1/c * diag(4))
	
	mmm = inv_Sigma %*% colSums(theta)
	
	return( mvrnorm(1, cov_mat %*% mmm, cov_mat))
}

# draw from full conditional of Sigma
update_Sigma = function(theta, mu, Psi, nu){
	
	n = dim(theta)[1]
	
	S = matrix(0, 4,4)
	
	for(s in 1:n){
		temp = theta[s, ]-mu
		S = S + temp %*% t(temp)
	}
	
	
	
	return(rinvwishart(nu + n, Psi + S))
}
