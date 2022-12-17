

file <- "https://www4.stat.ncsu.edu/~bjreich/ST740/HCDN_annual_max.RData"
load(url(file))

# remove NA's
good_rows = rowSums(is.na( Y ))
Y = Y[ which ( good_rows == 0) , ]

start_year = 1950
X = (seq(start_year, start_year + nyears-1, 1) - 1985) / 10
W = matrix(c(rep(1, nyears), X), nyears, 2)


max_iteration = 10000
burnin = 5000

seed = 1978
set.seed(seed)

source("mcmc.R")




max_iteration=100000
burnin = 25000


######################
out = MCMC(Y,
					W,
					# mcmc parameters
					max_iteration = max_iteration,
					burnin = burnin,
					# initial parameters with defaults
					mu0 = rep(0,4),
					Sigma0 = 1000*diag(4),
					# Sigma0 = diag(c(10000,10000,10,0.5)^2),
					# hyper parameters with defaults
					c = 100,
					Psi = diag(4),
					nu = 3)

out_noblock$acceptance_rate


ESS = rep(0,4)

for(loc in 1:n){
	samps <- mcmc(t(out$theta_trace[loc,,]),start=burnin,end=max_iteration,thin=1)
	
	ESS = ESS + effectiveSize(samps)
}

ESS = ESS / n
ESS

num_iter = max_iteration - burnin + 1
theta_df = data.frame(theta2 = c(out$theta_trace[1,2,burnin:max_iteration],
																 out$theta_trace[100,2,burnin:max_iteration],
																 out$theta_trace[200,2,burnin:max_iteration]),
											theta3 = c(out$theta_trace[1,3,burnin:max_iteration],
																 out$theta_trace[100,3,burnin:max_iteration],
																 out$theta_trace[200,3,burnin:max_iteration]),
											theta4 = c(out$theta_trace[1,4,burnin:max_iteration],
																 out$theta_trace[100,4,burnin:max_iteration],
																 out$theta_trace[200,4,burnin:max_iteration]),
											Location = c(rep("1",num_iter ),
																	 rep("100", num_iter ),
																	 rep("200", num_iter)),
											iteration=c(seq(burnin:max_iteration),
																	seq(burnin:max_iteration),
																	seq(burnin:max_iteration)))

theta2_plot = ggplot(theta_df, aes(y=theta2, x = iteration, col=Location))+ 
	geom_line(alpha = 0.5) +
	guides(colour = guide_legend(override.aes = list(size=4))) +
	ylab(TeX("$\\theta_{2}$")) +
	xlab("Iteration")+
	ggtitle(TeX("Trace of $\\theta_{2}$ at different locations"))+
	theme(text = element_text(size = 20))  
theta2_plot
ggsave("img/theta2_plot.png", theta2_plot)


theta3_plot = ggplot(theta_df, aes(y=theta3, x = iteration, col=Location))+ 
	geom_line(alpha = 0.5) +
	guides(colour = guide_legend(override.aes = list(size=4))) +
	ylab(TeX("$\\theta_{3}$")) +
	xlab("Iteration")+
	ggtitle(TeX("Trace of $\\theta_{3}$ at different locations"))+
	theme(text = element_text(size = 20))  
theta3_plot
ggsave("img/theta3_plot.png", theta3_plot)

theta4_plot = ggplot(theta_df, aes(y=theta4, x = iteration, col=Location))+ 
	geom_line(alpha = 0.5) +
	guides(colour = guide_legend(override.aes = list(size=4))) +
	ylab(TeX("$\\theta_{4}$")) +
	xlab("Iteration")+
	ggtitle(TeX("Trace of $\\theta_{4}$ at different locations"))+
	theme(text = element_text(size = 20))  
theta4_plot
ggsave("img/theta4_plot.png", theta4_plot)

library(maps)
library(mapproj)
library(viridis)
# s = s[ which ( good_rows == 0) , ]
map("state")

# points(s)
points(s[c(1, 100, 200),], pch= 18, col=2, cex=10)



make_map = function(Z, s, title, legend_title, file_name){
	dat <- data.frame(long=s[,1],lat=s[,2],Z=Z)
	map_plot = ggplot(dat, aes(long, lat)) +
		borders("state") +
		geom_point(aes(colour = Z), size = 4) +
		scale_colour_gradientn(legend_title, colours = viridis(10)) +
		coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
		xlab("")+ylab("")+labs(title=title) +theme(text = element_text(size = 20))
	
	ggsave(paste0("img/", file_name), map_plot)
	
	
}

posterior_means = rowMeans(out$theta_trace[,2,burnin:max_iteration])
make_map(posterior_means, s, TeX("Posterior mean of $\\theta_{2}$"), TeX("$\\theta_{2}$ mean"),"map_posterior_mean.png")

# posterior_sds = apply(out2$theta_trace[,2,burnin:max_iteration], c(1,2), sd)
# make_map(posterior_sds[,2], s, TeX("Posterior standard deviations of $\\theta_{2}$"), TeX("$\\theta_{2}$ sd"), "map_posterior_sds.png")

posterior_probs = rowMeans(out$theta_trace[,2,burnin:max_iteration] > 0)
make_map(posterior_probs, s, TeX("Posterior probability of $\\theta_{2}>0$"), TeX("$P(\\theta_{2}>0)$"),"map_posterior_prob.png")

summary(posterior_probs)
mean(posterior_probs > 0.5)




###############
# Plotting
###############
library(ggplot2)
library(latex2exp)





###
# Test different hyperparameters of mu and Sigma
###


out1 = MCMC(Y,
						W,
						# mcmc parameters
						max_iteration = max_iteration,
						burnin = burnin,
						# initial parameters with defaults
						mu0 = rep(0,4),
						Sigma0 = 1000*diag(4),
						# Sigma0 = diag(c(10000,10000,10,0.5)^2),
						# hyper parameters with defaults
						c = 100,
						Psi = diag(4),
						nu = 3)


out2 = MCMC(Y,
						W,
						# mcmc parameters
						max_iteration = max_iteration,
						burnin = burnin,
						# initial parameters with defaults
						mu0 = rep(0,4),
						Sigma0 = 1000*diag(4),
						# Sigma0 = diag(c(10000,10000,10,0.5)^2),
						# hyper parameters with defaults
						c = 1,
						Psi = diag(4),
						nu = 3)


out3 = MCMC(Y,
						W,
						# mcmc parameters
						max_iteration = max_iteration,
						burnin = burnin,
						# initial parameters with defaults
						mu0 = rep(0,4),
						Sigma0 = 1000*diag(4),
						# Sigma0 = diag(c(10000,10000,10,0.5)^2),
						# hyper parameters with defaults
						c = 1,
						Psi = 4 * diag(4),
						nu = 9)

num_iter = max_iteration - burnin + 1

mu_Sigma_df = data.frame(mu1 = c(out1$mu_trace[,1,burnin:max_iteration],
																	out2$mu_trace[,1,burnin:max_iteration],
																	out3$mu_trace[,1,burnin:max_iteration]),
												 Sigma11 = c(out1$Sigma_trace[1,1,burnin:max_iteration],
												 						out2$Sigma_trace[1,1,burnin:max_iteration],
												 						out3$Sigma_trace[1,1,burnin:max_iteration]),
												 Parameters = c(rep("(1)",num_iter ),
												 					 rep("(2)", num_iter ),
												 					 rep("(3)", num_iter)),
												 iteration=c(seq(burnin:max_iteration),
												 						seq(burnin:max_iteration),
												 						seq(burnin:max_iteration)))

mu_plot = ggplot(mu_Sigma_df, aes(y=mu1, x = iteration, col=Parameters))+ 
	geom_line(alpha = 0.5) +
	guides(colour = guide_legend(override.aes = list(size=4))) +
	ylab(TeX("$\\mu_{1}$")) +
	xlab("Iteration")+
	ggtitle(TeX("Trace of $\\mu_{1}$ under different prior parameters"))+
	theme(text = element_text(size = 20))  
mu_plot
ggsave("img/mu_plot.png", mu_plot)

Sigma_plot = ggplot(mu_Sigma_df, aes(y=Sigma11, x = iteration, col=Parameters))+ 
	geom_line(alpha = 0.5) +
	guides(colour = guide_legend(override.aes = list(size=4))) +
	ylab(TeX("$\\Sigma_{1,1}$")) +
	xlab("Iteration") +
	ggtitle(TeX("Trace of $\\Sigma_{1,1}$ under different prior parameters"))+
	theme(text = element_text(size = 20))  
Sigma_plot
ggsave("img/sigma_plot.png", Sigma_plot)
