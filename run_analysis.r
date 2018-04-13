#########################################################
## STAT 619 Term Project - Analysis scripts
## Andy Pohl - 30058005
## University of Calgary
## March 2017
#########################################################

##### Peliminaries #####
library(rjags)
library(R2WinBUGS)
library(coda)
library(BRugs)
library(latex2exp)

options(scipen=999) # set scientific notation for plotting clarity

##### Set workin directory #####
setwd("C:/Users/AndyPohl/Dropbox/PhD - Calgary/Course_work/STAT619/Project/Analysis/EFM")

### Source key functions and data ####

source("./code/get_data.r")
source("./code/generate_model_files.r")
source("./code/Bayes_support_functions.r")


#### Set save parameter ####
# if true then MCMC chain will be saved in the MCMC folder
saveMCMC = TRUE 

# if you have previously run the various models and saved the results you can load them as follows
load('./MCMC/fe_sim.RData')
load('./MCMC/re_sim_1.RData');load('./MCMC/re_sim_2.RData');
load('./MCMC/re_sim_3.RData');load('./MCMC/re_sim_4.RData');
load('./MCMC/re_sim_5.RData'); 

load('./MCMC/obs_sim.RData'); load('./MCMC/equivalent_prior_sim.RData');
load('./MCMC/naive_prior_sim.Rdata'); load('./MCMC/skeptical_prior_sim.RData'); 


load('./MCMC/gs_sim.RData'); load('./MCMC/gs_sim_2.RData')
load(file = './MCMC/regression_sim.RData')

all_re_models = list(re_1 = re_sim_1, re_2 = re_sim_2, re_3 = re_sim_3, re_4 = re_sim_4, re_5 = re_sim_5)

# if not you may run all models by sourcing the following:

source(file = "./code/run_models.r")

###################################################################################
### Generate results
###################################################################################

# The following code generates the various summary figures found within the manuscript.  
# Note that posterior quantites of interest can be examined via accessing the appropiate
# column of the $sims.matrix for that parameter.  Not every posterior quantity mentioned
# within the paper is accessed in the following code...

###################################################################################
### Fixed effects model ###
###################################################################################

# Assess convergence
# Plot chain...
par(mfrow = c(2,1))
plot((burnin+1):n.iter, fe_sim$sims.array[,1,1], type = 'l', col ='blue',
     main = TeX('Fixed Effects Model: Trace Plot for $\\theta$'),
     xlab = 'Index', 
     ylab = 'Value',
     ylim = c(-3,3))
lines((burnin+1):n.iter, fe_sim$sims.array[,2,1], type = 'l', col ='red')
lines((burnin+1):n.iter, fe_sim$sims.array[,3,1], type = 'l', col ='green')

gelman.plot(as.mcmc.list(lapply(as.data.frame(fe_sim$sims.array[,,1]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\theta$"), auto.layout = F)

head(fe_sim$sims.matrix)

# model checks
fe_sim$summary

# Posterior Density Plot
par(mfrow = c(1,1))
plotPost(fe_sim$sims.matrix[,1], showCurve = T, cenTend = "mean", disbCRL = T,dispHDI = F,
         xlab = TeX("$\\theta$"), main = TeX("Posterior Density for $\\theta$"))


###################################################################################
### Random Effects Model ###
###################################################################################

## choose a model index 1-5 to select the appropiate model in all_re_models list...
mdl_idx = 1

## Convergence Checks
    #trace plots
par(mfrow = c(2,1))
plot((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,1,4], type = 'l', col ='blue',
     main = TeX('Trace Plot for $\\mu$'),
     xlab = 'Index', 
     ylab = 'Value')
lines((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,2,4], type = 'l', col ='red')
lines((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,3,4], type = 'l', col ='green')

plot((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,1,2], type = 'l', col ='blue',
     main = TeX('Trace Plot for $\\tau^2$'),
     xlab = 'Index', 
     ylab = 'Value')
lines((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,2,2], type = 'l', col ='red')
lines((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,3,2], type = 'l', col ='green')

gelman.plot(as.mcmc.list(lapply(as.data.frame(all_re_models[[mdl_idx]]$sims.array[,,4]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\mu$"), auto.layout = F)

gelman.plot(as.mcmc.list(lapply(as.data.frame(all_re_models[[mdl_idx]]$sims.array[,,2]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\tau^2$"), auto.layout = F)
            

par(mfrow = c(3,3))
for(j in 1:9){
  plot((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,1,4+j], type = 'l', col ='blue',
       main = TeX(sprintf('Trace Plot for $\\theta_{%i}$',j)),
       xlab = 'Index', 
       ylab = 'Value')
  lines((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,2,4+j], type = 'l', col ='red')
  lines((burnin+1):n.iter, all_re_models[[mdl_idx]]$sims.array[,3,4+j], type = 'l', col ='green')
}

par(mfrow = c(3,3))
for(j in 1:9){
  gelman.plot(as.mcmc.list(lapply(as.data.frame(all_re_models[[mdl_idx]]$sims.array[,,4+j]), mcmc)),
              main = TeX(sprintf("Gelman-Rubian Diagonistic for $\\theta_{%i}$", j)), auto.layout = F)
}


## Sensitivity analysis
# the matrix sensitivity analysis contains the various summary statistics for the 5 
# random effects models investigated

sensitivity_analysis = matrix(NA, 15,5)
for(i in 1:5){
  #mu
  sensitivity_analysis[1,i] = mean(all_re_models[[i]]$sims.matrix[,4])
  sensitivity_analysis[2,i] = sd(all_re_models[[i]]$sims.matrix[,4])
  sensitivity_analysis[3,i] = quantile(all_re_models[[i]]$sims.matrix[,4], probs = c(0.025))
  sensitivity_analysis[4,i] = quantile(all_re_models[[i]]$sims.matrix[,4], probs = c(0.975))
  sensitivity_analysis[5,i] = mean(all_re_models[[i]]$sims.matrix[,4]<0)
  
  #tau2
  sensitivity_analysis[6,i] = mean(all_re_models[[i]]$sims.matrix[,2])
  sensitivity_analysis[7,i] = sd(all_re_models[[i]]$sims.matrix[,2])
  sensitivity_analysis[8,i] = quantile(all_re_models[[i]]$sims.matrix[,2], probs = c(0.025))
  sensitivity_analysis[9,i] = quantile(all_re_models[[i]]$sims.matrix[,2], probs = c(0.975))
  
  #dic
  sensitivity_analysis[10,i] = mean(all_re_models[[i]]$sims.matrix[,28])
  #model checks - p.values
  sensitivity_analysis[11,i] = mean(all_re_models[[i]]$sims.matrix[,25]) #max
  sensitivity_analysis[12,i] = mean(all_re_models[[i]]$sims.matrix[,26]) #min
  sensitivity_analysis[13,i] = mean(all_re_models[[i]]$sims.matrix[,27]) #range
  sensitivity_analysis[14,i] = mean(all_re_models[[i]]$sims.matrix[,14]) #mean
  sensitivity_analysis[15,i] = mean(all_re_models[[i]]$sims.matrix[,15]) #sd
}

### Create Forest Plot ###
nstudies = nrow(efm_data[efm_data$study_type==1,])
mdlIdx = 1
theta_samps = all_re_models[[mdlIdx]]$sims.list$theta # get thetas for the model 

openGraph( height=10 , width=10 )
par( mar=c(3.0,4.0,0.5,0.5) , mgp=c(2.0,0.7,0) )
muLim = quantile( theta_samps, 
                  probs=c(0.001,0.999) )
xLim = c( muLim[1] , muLim[2] + 0.5*(muLim[2]-muLim[1]) )
plot( -1,-1, 
      xlim=xLim ,
      xlab="Risk Difference (per 1000 people)" ,
      ylim=c(1,nstudies+2.5) , ylab="Study ID" ,
      yaxt="n" )
axis( side=2 , at=c(3:(nstudies+2),1) ,
      labels=c(as.character(paste0('R',1:nstudies)),"Overall") , las=1 )
abline( v=0 , lty="dotted" )
for ( sIdx in 1:nstudies ) {
  origInt = c(efm_data$RD[sIdx] - 1.96*sqrt(efm_data$var_RD[sIdx]), efm_data$RD[sIdx] + 1.96*sqrt(efm_data$var_RD[sIdx]))
  print(origInt)
  plotViolin( mcmcSamp=theta_samps[,sIdx] , plotHt=sIdx+2 ,
              cenTend = 'mean',
              dataVal=(efm_data$RD[sIdx]), 
              dataN=(efm_data$n[sIdx]) ,
              origInt = origInt,
              topOnly=TRUE, axisLim = xLim )
}
# overall
plotViolin(mcmcSamp = all_re_models[[mdlIdx]]$sims.list$mu, plotHt = 1, 
           topOnly = TRUE, cenTend = 'mean')
abline(h=2)


saveGraph(file = './figures/re_forest_plot', type = 'pdf')
saveGraph(file = './figures/re_forest_plot', type = 'png')

#######################################################################
### Observational Evidence - independent analysis
#######################################################################

## convergence checks
#trace plots
par(mfrow = c(2,1))
plot((burnin+1):n.iter, obs_sim$sims.array[,1,4], type = 'l', col ='blue',
     main = TeX('Trace Plot for $\\mu$'),
     xlab = 'Index', 
     ylab = 'Value')
lines((burnin+1):n.iter, obs_sim$sims.array[,2,4], type = 'l', col ='red')
lines((burnin+1):n.iter, obs_sim$sims.array[,3,4], type = 'l', col ='green')

plot((burnin+1):n.iter, obs_sim$sims.array[,1,2], type = 'l', col ='blue',
     main = TeX('Trace Plot for $\\tau^2$'),
     xlab = 'Index', 
     ylab = 'Value')
lines((burnin+1):n.iter, obs_sim$sims.array[,2,2], type = 'l', col ='red')
lines((burnin+1):n.iter, obs_sim$sims.array[,3,2], type = 'l', col ='green')


gelman.plot(as.mcmc.list(lapply(as.data.frame(obs_sim$sims.array[,,4]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\mu$"), auto.layout = F)

gelman.plot(as.mcmc.list(lapply(as.data.frame(obs_sim$sims.array[,,2]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\tau^2$"), auto.layout = F)

par(mfrow = c(3,3))
for(j in 1:17){
  plot((burnin+1):n.iter, obs_sim$sims.array[,1,4+j], type = 'l', col ='blue',
       main = TeX(sprintf('Trace Plot for $\\theta_{%i}$',j)),
       xlab = 'Index', 
       ylab = 'Value')
  lines((burnin+1):n.iter, obs_sim$sims.array[,2,4+j], type = 'l', col ='red')
  lines((burnin+1):n.iter, obs_sim$sims.array[,3,4+j], type = 'l', col ='green')
}

par(mfrow = c(3,3))
for(j in 1:17){
  gelman.plot(as.mcmc.list(lapply(as.data.frame(obs_sim$sims.array[,,4+j]), mcmc)),
              main = TeX(sprintf("Gelman-Rubian Diagonistic for $\\theta_{%i}$", j)), auto.layout = F)
}

# create forest plots
nstudies = nrow(efm_data[efm_data$study_type!=1,])
mdlIdx = 1
theta_samps = obs_sim$sims.list$theta # get thetas for the model 

openGraph( height=10 , width=10 )
par( mar=c(3.0,4.0,0.5,0.5) , mgp=c(2.0,0.7,0) )
muLim = quantile( theta_samps, 
                  probs=c(0.001,0.999) )
xLim = c( muLim[1] , muLim[2] + 0.5*(muLim[2]-muLim[1]) )
plot( -1,-1, 
      xlim=xLim ,
      xlab="Risk Difference (per 1000 people)" ,
      ylim=c(1,nstudies+2.5) , ylab="Study ID" ,
      yaxt="n" )
axis( side=2 , at=c(3:(nstudies+2),1) ,
      labels=c(as.character(paste0('C',1:7)),as.character(paste0('B',1:10)),"Overall") , las=1 )
abline( v=0 , lty="dotted" )
for ( sIdx in 1:nstudies ) {
  origInt = c(efm_data$RD[sIdx+9] - 1.96*sqrt(efm_data$var_RD[sIdx+9]), efm_data$RD[sIdx+9] + 1.96*sqrt(efm_data$var_RD[sIdx+9]))
  print(origInt)
  plotViolin( mcmcSamp=theta_samps[,sIdx] , plotHt=sIdx+2 ,
              cenTend = 'mean',
              dataVal=(efm_data$RD[sIdx+9]), 
              dataN=(efm_data$n[sIdx+9]) ,
              origInt = origInt,
              topOnly=TRUE, axisLim = xLim )
}
# overall
plotViolin(mcmcSamp = obs_sim$sims.list$mu, plotHt = 1, 
           topOnly = TRUE, cenTend = 'mean')
abline(h=2)
saveGraph(file = './figures/obs_forest_plot', type = 'pdf')
saveGraph(file = './figures/obs_forest_plot', type = 'png')


#######################################################################
### Observational Evidence - Influence of prior
#######################################################################

##### make prior comparision plot

par(mfrow = c(1,4))
x = seq(-6,3,0.001)
likelihood = dnorm(x, mean = -0.25, sd = 1.12)
ref_prior = dnorm(x, mean = 0, sd = sqrt(10^6))
ref_posterior = dnorm(x, mean = mean(all_re_models[[1]]$sims.list$mu), 
                      sd = sd(all_re_models[[1]]$sims.list$mu))
plot(x, likelihood, type = 'l', ylim = c(0,1), xlim = c(-6,3), lty = 4,
     xlab = 'RD (per 1000 people)', ylab = 'Density', main = 'Vague Prior')
lines(x, ref_prior, type = 'l', lty = 2, col = 'darkgrey')
lines(x, ref_posterior, lty = 1, col = 'blue')
abline(v = 0, lty = 3)

skeptical_prior_dens = dnorm(x, mean = -1.637937,sd=sqrt(5.062186))
skeptical_posterior_dens = dnorm(x, mean = mean(skeptical_prior_sim$sims.list$mu), 
                                 sd = sd(skeptical_prior_sim$sims.list$mu))
plot(x, likelihood, type = 'l', ylim = c(0,1), xlim = c(-6,3), lty = 4,
     xlab = 'RD (per 1000 people)', ylab = 'Density', col = 'red', main = 'Skeptical Prior')
lines(x, skeptical_prior_dens, type = 'l', lty = 2, col = 'darkgrey')
lines(x, skeptical_posterior_dens, lty = 1, col = 'blue')
abline(v = 0, lty = 3)

equivalent_prior_dens = dnorm(x, mean = -1.637937,sd=sqrt(1.265547))
equivalent_posterior_dens = dnorm(x, mean = mean(equivalent_prior_sim$sims.list$mu), 
                                  sd = sd(equivalent_prior_sim$sims.list$mu))
plot(x, likelihood, type = 'l', ylim = c(0,1), xlim = c(-6,3), lty = 4,
     xlab = 'RD (per 1000 people)', ylab = 'Density', col = 'red', main = 'Equivalent Prior')
lines(x, equivalent_prior_dens, type = 'l', lty = 2, col = 'darkgrey')
lines(x, equivalent_posterior_dens, lty = 1, col = 'blue')
abline(v = 0, lty = 3)

naive_prior_dens = dnorm(x, mean = -1.637937,sd=sqrt(0.2055758))
naive_posterior_dens = dnorm(x, mean = mean(naive_prior_sim$sims.list$mu), 
                             sd = sd(naive_prior_sim$sims.list$mu))
plot(x, likelihood, type = 'l', ylim = c(0,1), xlim = c(-6,3), lty = 4,
     xlab = 'RD (per 1000 people)', ylab = 'Density', col = 'red', main = 'Naive Prior')
lines(x, naive_prior_dens, type = 'l', lty = 2, col = 'darkgrey')
lines(x, naive_posterior_dens, lty = 1, col = 'blue')
abline(v = 0, lty = 3)


### get results ##
model_list = list(all_re_models[[1]], skeptical_prior_sim, equivalent_prior_sim, naive_prior_sim)
sensitivity_analysis = matrix(NA, 10,length(model_list))
for(i in 1:length(model_list)){
  #mu
  sensitivity_analysis[1,i] = mean(model_list[[i]]$sims.list$mu)
  sensitivity_analysis[2,i] = sd(model_list[[i]]$sims.list$mu)
  sensitivity_analysis[3,i] = quantile(model_list[[i]]$sims.list$mu, probs = c(0.025))
  sensitivity_analysis[4,i] = quantile(model_list[[i]]$sims.list$mu, probs = c(0.975))
  sensitivity_analysis[5,i] = mean(model_list[[i]]$sims.list$mu<0)
  
  #tau2
  sensitivity_analysis[6,i] = mean(model_list[[i]]$sims.list$var)
  sensitivity_analysis[7,i] = sd(model_list[[i]]$sims.list$var)
  sensitivity_analysis[8,i] = quantile(model_list[[i]]$sims.list$var, probs = c(0.025))
  sensitivity_analysis[9,i] = quantile(model_list[[i]]$sims.list$var, probs = c(0.975))
  
  sensitivity_analysis[10, i] = mean(model_list[[i]]$sims.list$deviance)
}

#######################################################################
### Generalised synthesis
#######################################################################

### Plot traceplots
## Mu
par(mfrow=c(2,1))
plot(y = gs_sim$sims.array[,1,1],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\mu$'))
lines(x=5001:50000, y=gs_sim$sims.array[,2,1], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim$sims.array[,3,1], col = 'green', type = 'l')

## tau
plot(y = gs_sim$sims.array[,1,35],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\tau^2$'))
lines(x=5001:50000, y=gs_sim$sims.array[,2,35], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim$sims.array[,3,35], col = 'green', type = 'l')

gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim$sims.array[,,1]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\mu$"), auto.layout = F)

gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim$sims.array[,,35]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\tau^2$"), auto.layout = F)

### Plot traceplots
par(mfrow=c(3,1))
## theta1
plot(y = gs_sim$sims.array[,1,3],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\theta_1$'))
lines(x=5001:50000, y=gs_sim$sims.array[,2,3], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim$sims.array[,3,3], col = 'green', type = 'l')
## theta2
plot(y = gs_sim$sims.array[,1,4],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\theta_2$'))
lines(x=5001:50000, y=gs_sim$sims.array[,2,4], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim$sims.array[,3,4], col = 'green', type = 'l')
## theta3
plot(y = gs_sim$sims.array[,1,5],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\theta_3$'))
lines(x=5001:50000, y=gs_sim$sims.array[,2,5], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim$sims.array[,3,5], col = 'green', type = 'l')

gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim$sims.array[,,3]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\theta_1$"), auto.layout = F)
gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim$sims.array[,,4]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\theta_2$"), auto.layout = F)
gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim$sims.array[,,5]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\theta_3$"), auto.layout = F)

### Create forest plot

nstudies = nrow(efm_data)
nstudytypes = length(unique(efm_data$study_type))
nrct = 9
ncc = 7
nba = 10

openGraph( height=40 , width=20 )
par( mar=c(3.0,5.5,0.5,0.5) , mgp=c(2.0,0.7,0) )
muLim = quantile( gs_sim$sims.list$psy, 
                  probs=c(0.001,0.999) )
xLim = c( muLim[1] , muLim[2] + 0.5*(muLim[2]-muLim[1]) )
plot( -1,-1, 
      xlim=xLim ,
      xlab="Risk Difference (per 1000 people)" , ylab = '',
      ylim=c(1,33.5),
      yaxt="n" )
axis( side=2 , at=c(1,3:13, 15:22, 24:33) ,
      labels=c("Overall", "Overall BA", as.character(paste0('B',nba:1)),
               "Overall CC", as.character(paste0('C',ncc:1)), 
               "Overall RCT", as.character(paste0('R',nrct:1))) , las=1 )
abline( v=0 , lty="dotted" )
abline(h=2)
abline(h = 14, lty = 'dashed', colour = 'darkgrey')
abline(h = 23, lty = 'dashed', colour = 'darkgrey')
plotViolin(mcmcSamp = gs_sim$sims.matrix[,1], plotHt = 1, 
           topOnly = TRUE, cenTend = 'mean')
plotViolin(mcmcSamp = gs_sim$sims.matrix[,5], plotHt = 3, 
           topOnly = TRUE, cenTend = 'mean')

plotViolin(mcmcSamp = gs_sim$sims.matrix[,4], plotHt = 15, 
           topOnly = TRUE, cenTend = 'mean')
plotViolin(mcmcSamp = gs_sim$sims.matrix[,3], plotHt = 24, 
           topOnly = TRUE, cenTend = 'mean')
for(sIdx in 1:nrct){
  
  ht = 33-sIdx+1 
  print(ht)
  print(sIdx)
  #origInt = c(efm_data$RD[sIdx] - 1.96*sqrt(efm_data$var_RD[sIdx]), efm_data$RD[sIdx] + 1.96*sqrt(efm_data$var_RD[sIdx]))
  #print(origInt)
  plotViolin( mcmcSamp=gs_sim$sims.matrix[,5+sIdx] , plotHt=ht ,
              cenTend = 'mean',
              dataVal=(efm_data$RD[sIdx]), 
              dataN=(efm_data$n[sIdx]) ,
              #origInt = origInt,
              topOnly=TRUE, axisLim = xLim )
}
for(sIdx in 1:ncc){
  
  ht = 22-sIdx+1 
  print(ht)
  print(sIdx)
  #origInt = c(efm_data$RD[sIdx] - 1.96*sqrt(efm_data$var_RD[sIdx]), efm_data$RD[sIdx] + 1.96*sqrt(efm_data$var_RD[sIdx]))
  #print(origInt)
  plotViolin( mcmcSamp=gs_sim$sims.matrix[,14+sIdx] , plotHt=ht ,
              cenTend = 'mean',
              dataVal=(efm_data$RD[sIdx+9]), 
              dataN=(efm_data$n[sIdx+9]) ,
              #origInt = origInt,
              topOnly=TRUE, axisLim = xLim )
}

for(sIdx in 1:nba){
  
  ht = 13-sIdx+1 
  print(ht)
  print(sIdx)
  #origInt = c(efm_data$RD[sIdx] - 1.96*sqrt(efm_data$var_RD[sIdx]), efm_data$RD[sIdx] + 1.96*sqrt(efm_data$var_RD[sIdx]))
  #print(origInt)
  plotViolin( mcmcSamp=gs_sim$sims.matrix[,21+sIdx] , plotHt=ht ,
              cenTend = 'mean',
              dataVal=(efm_data$RD[sIdx+16]), 
              dataN=(efm_data$n[sIdx+16]) ,
              #origInt = origInt,
              topOnly=TRUE, axisLim = xLim )
}
saveGraph(file = './figures/gs_forest_plot', type = 'pdf')
saveGraph(file = './figures/gs_forest_plot', type = 'png')


## Generalised synthesis model 2 with half cauchy priors

### Plot traceplots
## Mu
par(mfrow=c(2,1))
plot(y = gs_sim_2$sims.array[,1,1],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\mu$'))
lines(x=5001:50000, y=gs_sim_2$sims.array[,2,1], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim_2$sims.array[,3,1], col = 'green', type = 'l')

## tau
plot(y = gs_sim_2$sims.array[,1,34],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\tau^2$'))
lines(x=5001:50000, y=gs_sim_2$sims.array[,2,34], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim_2$sims.array[,3,34], col = 'green', type = 'l')

gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim_2$sims.array[,,1]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\mu$"), auto.layout = F)

gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim_2$sims.array[,,34]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\tau^2$"), auto.layout = F)

### Plot traceplots
par(mfrow=c(3,1))
## theta1
plot(y = gs_sim_2$sims.array[,1,2],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\theta_1$'))
lines(x=5001:50000, y=gs_sim_2$sims.array[,2,2], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim_2$sims.array[,3,2], col = 'green', type = 'l')
## theta2
plot(y = gs_sim_2$sims.array[,1,3],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\theta_2$'))
lines(x=5001:50000, y=gs_sim_2$sims.array[,2,3], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim_2$sims.array[,3,3], col = 'green', type = 'l')
## theta3
plot(y = gs_sim_2$sims.array[,1,4],x=5001:50000, col = 'red', type = 'l',
     ylab = 'Value', xlab = 'index', main = TeX('Trace Plot for $\\theta_3$'))
lines(x=5001:50000, y=gs_sim_2$sims.array[,2,4], col = 'blue', type = 'l')
lines(x=5001:50000, y=gs_sim_2$sims.array[,3,4], col = 'green', type = 'l')

gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim_2$sims.array[,,2]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\theta_1$"), auto.layout = F)
gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim_2$sims.array[,,3]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\theta_2$"), auto.layout = F)
gelman.plot(as.mcmc.list(lapply(as.data.frame(gs_sim_2$sims.array[,,4]), mcmc)),
            main = TeX("Gelman-Rubian Diagonistic for $\\theta_3$"), auto.layout = F)



