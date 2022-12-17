library(ggplot2)
library(maps)
library(mapproj)
library(viridis)


source("hw3/update_params.R")

file <- "https://www4.stat.ncsu.edu/~bjreich/ST740/HCDN_annual_max.RData"
load(url(file))

# remove NA's
good_rows = rowSums(is.na( Y ))
Y = Y[ which ( good_rows == 0) , ]

Z = log(Y+1)
start_year = 1950
X = (seq(start_year, start_year + nyears-1, 1) - 1985) / 10
W = matrix(c(rep(1, nyears), X), nyears, 2)

dat <- data.frame(long=s[,1],lat=s[,2],Z=Z1)
ggplot(dat, aes(long, lat)) +
	borders("state") +
	geom_point(aes(colour = Z)) +
	scale_colour_gradientn(colours = viridis(10)) +
	coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
	xlab("")+ylab("")+labs(title="log streamflow in 2021")



# sample from multivariate normal
library(MASS)

seed = 1978
set.seed(seed)
max_iteration = 10000
burnin = 5000


n = dim(Z)[1]
m = dim(Z)[2]


###
# initial parameters
###
mu = c(0,0)
Sigma = diag(2)
sigmasq = 1
beta = mvrnorm(n, mu, Sigma)

# hyper parameters
c = 100
a = 0.5
b = 0.5
Psi = diag(2)
nu = 3

# start trace matrices
mu_trace = array(0, c(1, 2, max_iteration))
mu_trace[,,1] = mu

Sigma_trace = array(0, c(2,2, max_iteration))
Sigma_trace[,,1] = Sigma

sigmasq_trace = array(0, c(1, max_iteration))
sigmasq_trace[,1] = sigmasq

beta_trace=array(0,c(n,2,max_iteration))
beta_trace[,,1] = beta

###
# run MCMC
###
for(iteration in 2:max_iteration){
	
	if(iteration %% 200 == 0){print(iteration)}
	
	# update parameters
	mu = update_mu(Sigma, beta, c)
	sigmasq = update_sigmasq(Z, W, beta, a, b)
	Sigma = update_Sigma(beta, mu, Psi, nu)
	
	for(s in 1:n){beta[s,] = update_beta(Z[s,], W, mu, Sigma, sigmasq)}
	
	# save current value to trace
	mu_trace[,,iteration] = mu
	Sigma_trace[,,iteration] = Sigma
	sigmasq_trace[,iteration] = sigmasq
	beta_trace[,,iteration] = beta
}

###
# trace plots
###
library(ggplot2)
library(latex2exp)


make_trace_plot = function(trace, ylab, file_name){
	
	trace_plot = ggplot(data.frame(y=trace), aes(x=1:length(trace), y=y)) +
			geom_line() +
			ylab(TeX(ylab)) +
			xlab("iteration") + 
			ggtitle(TeX(paste0("Trace plot for ", ylab))) +
			theme(text = element_text(size = 20))  
	
	ggsave(paste0("hw3/img/", file_name), trace_plot)
}

make_trace_plot(mu_trace[1,1,burnin:max_iteration], "$\\mu_{1}$", "trace_mu1.png")
make_trace_plot(mu_trace[1,2,burnin:max_iteration], "$\\mu_{2}$", "trace_mu2.png")

make_trace_plot(Sigma_trace[1,1,burnin:max_iteration], "$\\Sigma_{1,1}$", "trace_Sigma11.png")
make_trace_plot(Sigma_trace[2,1,burnin:max_iteration], "$\\Sigma_{2,1}$", "trace_Sigma21.png")
make_trace_plot(Sigma_trace[1,2,burnin:max_iteration], "$\\Sigma_{1,2}$", "trace_Sigma12.png")
make_trace_plot(Sigma_trace[2,2,burnin:max_iteration], "$\\Sigma_{2,2}$", "trace_Sigma22.png")

make_trace_plot(sigmasq_trace[1,burnin:max_iteration], "$\\sigma^{2}$", "trace_sigmasq.png")

make_trace_plot(beta_trace[1,1,burnin:max_iteration], "$\\beta_{0,1}$", "trace_beta01.png")
make_trace_plot(beta_trace[1,2,burnin:max_iteration], "$\\beta_{1,1}$", "trace_beta11.png")




###
# map plots
###

file <- "https://www4.stat.ncsu.edu/~bjreich/ST740/HCDN_annual_max.RData"
load(url(file))

make_map = function(Z, s, title, legend_title, file_name){
	dat <- data.frame(long=s_map[,1],lat=s_map[,2],Z=Z)
	map_plot = ggplot(dat, aes(long, lat)) +
		borders("state") +
		geom_point(aes(colour = Z), size = 4) +
		scale_colour_gradientn(legend_title, colours = viridis(10)) +
		coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
		xlab("")+ylab("")+labs(title=title) 
	
	ggsave(paste0("hw3/img/", file_name), map_plot)
	
	
}

posterior_means = rowMeans(beta_trace[,2,burnin:max_iteration])

s_map = s[ which ( good_rows == 0) , ]
make_map(posterior_means, s_map, TeX("Posterior mean of $\\beta_{1}$"), TeX("$\\beta_{1}$ mean"),"map_posterior_mean.png")

posterior_sds = apply(beta_trace[,burnin:max_iteration,], c(1,2), sd)
make_map(posterior_sds[,2], s_map, TeX("Posterior standard deviations of $\\beta_{1}$"), TeX("$\\beta_{1}$ sd"), "map_posterior_sds.png")

posterior_probs = rowMeans(beta_trace[,2,burnin:max_iteration] > 0)
make_map(posterior_probs, s_map, TeX("Posterior probability of $\\beta_{1}>0$"), TeX("$P(\\beta_{1}>0)$"),"map_posterior_prob.png")

summary(posterior_means)
