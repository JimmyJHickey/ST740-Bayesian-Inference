library(microbiome)
library(latex2exp)


data(dietswap)
Y_all <- otu_table(dietswap) # all data
X_all <- sample_data(dietswap)
Y <- Y_all[,X_all[,7]==1] # Extract baseline data
X <- X_all[X_all[,7]==1,3]

Y_no = Y[ rowSums(Y) != 0, ]
Y = Y_no
source("mcmc.R")
source("convergence_metrics.R")



N_vec = colSums(Y)

props = array(0, dim = c(m, n))

for(i in 1:n){
	props[,i] = Y[,i] / N_vec[i]
}

max_iteration = 2000
burnin = 200

mstar = mean(rowSds(props)^2)
vstar = sd(rowSds(props)^2)^2
a = (mstar^2 + 2 * vstar) / vstar
b = ( mstar^3 + mstar * vstar) / vstar

source("mcmc.R")

model1 = MCMC(Y,
					max_iteration = max_iteration,
					burnin = burnin,
					# if M_fixed is TRUE then M will be set to M0
					# otherwise its distribution will be estimated
					est_M = TRUE,
					MH = 0.2,
					# test hyper parameters
					# m0 = 50,
					# ssq = 1,
					# a = 1,
					# b = 1,
					m0 = mean(rowSums(props)),
					ssq = sd(rowSums(props))^2,
					a = a,
					b = b,
					seed = 1978)

model1_ESS = calc_ESS(model1$theta_trace, max_iteration, burnin)
mean(summary(model1_ESS)[4,])
model1$WAIC
# model1_geweke = calc_geweke(model1$theta_trace, max_iteration, burnin)

pdf(file="img/trace_plot_Aeromans.pdf")
plot(model1$theta_trace[1,2,200:2000], type = 'l', col="blue", xlab="iteration", ylab="Aeromans probability trace", cex.axis=1.25, cex.lab=1.5)
lines(model1$theta_trace[3,2,200:2000], type = 'l', xlab="iteration", ylab="Aeromans trace", col="orange")
lines(model1$theta_trace[2,2,200:2000], type = 'l', xlab="iteration", ylab="Aeromans trace", col="red")
dev.off()


pdf(file="img/trace_plot_Gemella.pdf")
plot(model1$theta_trace[1,61,200:2000], type = 'l', col="blue", xlab="iteration", ylab="Gemella probability trace", cex.axis=1.25, cex.lab=1.5)
lines(model1$theta_trace[3,61,200:2000], type = 'l', xlab="iteration", ylab="Gemella trace", col="orange")
lines(model1$theta_trace[2,61,200:2000], type = 'l', xlab="iteration", ylab="Gemella trace", col="red")
dev.off()


model2 = MCMC(Y,
						 max_iteration = max_iteration,
						 burnin = burnin,
						 # if M_fixed is TRUE then M will be set to M0
						 # otherwise its distribution will be estimated
						 est_M = FALSE,
						 M0 = 0,
						 seed = 1978)
model2$WAIC



model2_ESS = calc_ESS(model2$theta_trace, max_iteration, burnin)

# model2_geweke = calc_geweke(model2$theta_trace, max_iteration, burnin)

# model 3
alpha = rowSums(Y) / sum(Y)

theta = matrix(rep(alpha, n), n , m)

model3_ll = log_like(Y, theta)
model3_ll2 = model3_ll^2

mn_logpdf  <- model3_ll
# 0 variance because theta_i is fixed
var_logpdf <- model3_ll2 - model3_ll^2
pW         <- sum(var_logpdf)
WAIC       <- list(WAIC=-2*sum(mn_logpdf)+2*pW,pW=pW)
WAIC$WAIC




####################
# posterior pred check for model 2
####################
S = 1000
theta_thin = model2$theta_trace[,,seq(burnin, max_iteration, by = (max_iteration - burnin) / (S - 1))]
Y_tilde = array(0, dim = c(n, m, S))

n = dim(Y)[2]
m = dim(Y)[1]

for(i in 1:n){
	for(s in 1:S){
		Y_tilde[i,,s] = rmultinom(1, N_vec[i], theta_thin[i,,s])
	}
}

interval_check = array(0, dim = c(n,m))

for(i in 1:n){
	for(j in 1:m){
		interval_check[i,j] = (  quantile(Y_tilde[i,j,], 0.025) <= Y[j,i] &&  Y[j,i] <= quantile(Y_tilde[i,j,], 0.975) )
	}
	
}


plot_post_pred = function(Y_tilde, Y, i, j){
	pdf(file=paste("img/ppc_", i, "_", j, ".pdf", sep=''))
	hist(Y_tilde[i,j,], breaks=31, main=print(paste("Posterior predictive of ", rownames(Y[j]) , "\nfor patient ", i)),
			 xlab="Posterior predictive samples", cex.axis=1.25, cex.lab=1.5)
	abline(v = Y[j,i], lwd = 4, col = "blue")
	abline(v = quantile(Y_tilde[i,j,], 0.025), lwd = 4, col="black", lty="dotted")
	abline(v = quantile(Y_tilde[i,j,], 0.975), lwd = 4, col="black", lty="dotted")
	dev.off()
}

plot_post_pred(Y_tilde, Y, 1, 20)
plot_post_pred(Y_tilde, Y, 10, 1)
plot_post_pred(Y_tilde, Y, 10, 100)
plot_post_pred(Y_tilde, Y, 3, 89)



#############
# check difference between nationalities
#############
AAM_ind = (X == "AAM")
n_AAM = sum(AAM_ind)
AFR_ind = (X == "AFR")
n_AFR = sum(AFR_ind)

delta = array(0, dim = c(m, max_iteration-burnin))
delta_interval = array(0, dim = c(m,2))
delta_interval_check = array(0, dim = c(m))

for(j in 1:m){
	for(iteration in burnin:(max_iteration-1)) {
		delta[j, iteration - burnin] = sum(model2$theta_trace[ AAM_ind , j, iteration]) / n_AAM - 
								sum(model2$theta_trace[ AFR_ind , j, iteration]) / n_AFR 
	}
	
	delta_interval[j,1] = quantile( delta[j,] , 0.025)
	delta_interval[j,2] = quantile( delta[j,] , 0.975)
	delta_interval_check[j] = ( delta_interval[j,1] <= 0 & 0 <= delta_interval[j,2] ) 
	
}
sum(delta_interval_check)
delta_interval_length = delta_interval[,2] - delta_interval[,1]

sum(delta_interval_check[rowSums(Y) > 300])
sort(rowSums(Y))

sum(rowSums(Y) > 1000)


pdf(file=paste("img/delta.pdf", sep=''))
plot(seq(1,m), rep(0,m), type='l', lwd = 5, ylim=c(-0.0025, 0.0025),
		 cex.axis = 1.3, cex.lab = 1.35,
		 xlab="Taxa Number", ylab=TeX(r'($\Delta_{i}$ interval)'), main="Taxa difference posterior intervals")
for(j in 1:m){
	color = if(delta_interval_check[j]) "red" else "blue"
	
	segments(x0 = j, x1 = j, y0 = delta_interval[j,1], y1 = delta_interval[j,2], col = color, lwd = 4)
}
dev.off()




