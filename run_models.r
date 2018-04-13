#########################################################
## STAT 619 Term Project - Run models...
## Andy Pohl - 30058005
## University of Calgary
## March 2017
#########################################################

# MCMC settings
n.chains = 3
burnin = 5000
n.iter = 50000

#########################################################
## Fixed effects model
#########################################################

fe_model = "
model{
# prior
theta ~ dnorm(0, 0.000001)

# likelihood
for(i in 1:N){
p[i] <- 1/(sigma2[i])
y[i] ~ dnorm(theta, p[i])
y.rep[i] <- replicate.post(y[i])
}

## Posterior Checks

y.min <- ranked( y[1:N], 1 ) 
y.max <- ranked( y[1:N], N ) 
y.rep.min <- ranked( y.rep[1:N], 1 ) 
y.rep.max <- ranked( y.rep[1:N], N ) 
y.rep.min.dif <- y.min - y.rep.min 
y.rep.max.dif <- y.rep.max - y.max 

pred.y.min <- step( y.rep.min.dif ) 
pred.y.max <- step( y.rep.max.dif ) 

# Examine range 
y.rep.range.dif <- ( y.rep.max - y.rep.min ) - ( y.max - y.min ) 
pred.y.range <- step( y.rep.range.dif ) 

# Examine individual observations 
for( i in 1:N ) { 
y.rep.dif[i] <- y.rep[i] - y[i] 
pred.y[i] <- step( y.rep.dif[i] )
} 

# Examine sample and replicate sample mean and standard deviation 

pred.mean <- step( mean( y.rep[1:N] ) - mean( y[1:N] ) ) 
pred.sd <- step( sd(y.rep[1:N] ) - sd( y[1:N] ) )

}

}"

cat(fe_model,file="./modelFiles/fixed_effects.txt")
fe_sim = bugs(
  data = list(N = nrow(efm_data[efm_data$study_type==1,]), y = efm_data$RD[efm_data$study_type==1],  sigma2 = efm_data$var_RD[efm_data$study_type==1]),
  inits = list(list(theta = -10), list(theta = 10), list(theta = 0)),
  model.file = "./modelFiles/efm_fixed_effects.txt",
  parameters.to.save = c("theta", "pred.mean", "pred.sd", "pred.y", "pred.y.max", "pred.y.min", "pred.y.range"),
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = burnin,
  n.thin =1,
  bugs.directory = 'c:/Program Files (x86)/OpenBUGS/OpenBUGS323/', # suplement with path to bugs directory
  working.directory = getwd(),
  program = 'OpenBUGS')


#########################################################
## Random Effects Models ##
#########################################################

# Mdl 1: IG(0.001, 0.001)
re_model_1 = "
model{
# prior
tau ~ dgamma(0.001, 0.001)
var <- 1/tau
sd <- sqrt(var)

mu ~ dnorm(0, 0.000001)


# likelihood
for(i in 1:N){
p[i] <- 1/(sigma2[i])
y[i] ~ dnorm(theta[i], p[i])
theta[i] ~ dnorm (mu, tau)
y.rep[i] <- replicate.post(y[i])
}

# Posterior Checks
y.min <- ranked( y[1:N], 1 ) 
y.max <- ranked( y[1:N], N ) 
y.rep.min <- ranked( y.rep[1:N], 1 ) 
y.rep.max <- ranked( y.rep[1:N], N ) 
y.rep.min.dif <- y.min - y.rep.min 
y.rep.max.dif <- y.rep.max - y.max 

pred.y.min <- step( y.rep.min.dif ) 
pred.y.max <- step( y.rep.max.dif ) 

# Examine range 
y.rep.range.dif <- ( y.rep.max - y.rep.min ) - ( y.max - y.min ) 
pred.y.range <- step( y.rep.range.dif ) 

# Examine individual observations 
for( i in 1:N ) { 
y.rep.dif[i] <- y.rep[i] - y[i] 
pred.y[i] <- step( y.rep.dif[i] )
} 

# Examine sample and replicate sample mean and standard deviation 
pred.mean <- step( mean( y.rep[1:N] ) - mean( y[1:N] ) ) 
pred.sd <- step( sd(y.rep[1:N] ) - sd( y[1:N] ) )
}"
cat(re_model_1,file="./modelFiles/efm_re_1.txt")
re_sim_1 = bugs(
  data = list(N = nrow(efm_data[efm_data$study_type==1,]), y = efm_data$RD[efm_data$study_type==1],  sigma2 = efm_data$var_RD[efm_data$study_type==1]),
  inits = list(list(tau = 0.2, mu = -10),
               list(tau = 10, mu = 0),
               list(tau =  20, mu = 10)),
  model.file = "C:/Users/AndyPohl/Dropbox/PhD - Calgary/Course_work/STAT619/Project/Analysis/EFM/modelFiles/efm_re_1.txt",
  parameters.to.save = c("tau", "var", "sd", "mu", "theta", "pred.mean", "pred.sd", "pred.y", "pred.y.max", "pred.y.min", "pred.y.range"),
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = burnin,
  n.thin =1,
  bugs.directory = 'c:/Program Files (x86)/OpenBUGS/OpenBUGS323/',
  working.directory = getwd(),
  program = 'OpenBUGS')
save(re_sim_1, file = './MCMC/re_sim_1.RData')

##### Mdl 2: IG(0.1, 0.1) #####
re_model_2 = "
model{
# prior
tau ~ dgamma(0.1, 0.1)
var <- 1/tau
sd <- sqrt(var)

mu ~ dnorm(0, 0.000001)


# likelihood
for(i in 1:N){
p[i] <- 1/(sigma2[i])
y[i] ~ dnorm(theta[i], p[i])
theta[i] ~ dnorm (mu, tau)
y.rep[i] <- replicate.post(y[i])
}

# Posterior Checks
y.min <- ranked( y[1:N], 1 ) 
y.max <- ranked( y[1:N], N ) 
y.rep.min <- ranked( y.rep[1:N], 1 ) 
y.rep.max <- ranked( y.rep[1:N], N ) 
y.rep.min.dif <- y.min - y.rep.min 
y.rep.max.dif <- y.rep.max - y.max 

pred.y.min <- step( y.rep.min.dif ) 
pred.y.max <- step( y.rep.max.dif ) 

# Examine range 
y.rep.range.dif <- ( y.rep.max - y.rep.min ) - ( y.max - y.min ) 
pred.y.range <- step( y.rep.range.dif ) 

# Examine individual observations 
for( i in 1:N ) { 
y.rep.dif[i] <- y.rep[i] - y[i] 
pred.y[i] <- step( y.rep.dif[i] )
} 

# Examine sample and replicate sample mean and standard deviation 
pred.mean <- step( mean( y.rep[1:N] ) - mean( y[1:N] ) ) 
pred.sd <- step( sd(y.rep[1:N] ) - sd( y[1:N] ) )
}"
cat(re_model_2,file="./modelFiles/efm_re_2.txt")
re_sim_2 = bugs(
  data = list(N = nrow(efm_data[efm_data$study_type==1,]), y = efm_data$RD[efm_data$study_type==1],  sigma2 = efm_data$var_RD[efm_data$study_type==1]),
  inits = list(list(tau = 0.2, mu = -10),
               list(tau = 10, mu = 0),
               list(tau =  20, mu = 10)),
  model.file = "C:/Users/AndyPohl/Dropbox/PhD - Calgary/Course_work/STAT619/Project/Analysis/EFM/modelFiles/efm_re_2.txt",
  parameters.to.save = c("tau", "var", "sd", "mu", "theta", "pred.mean", "pred.sd", "pred.y", "pred.y.max", "pred.y.min", "pred.y.range"),
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = burnin,
  n.thin =1,
  bugs.directory = 'c:/Program Files (x86)/OpenBUGS/OpenBUGS323/',
  working.directory = getwd(),
  program = 'OpenBUGS')
save(re_sim_2, file = './MCMC/re_sim_2.RData')


##### Mdl 3: Half Normal #####
re_model_3 = "
## Half normal prior
model{
# prior
sd ~ dnorm(0,0.001)I(0,)
var <- pow(sd, 2)
tau <- 1/var

mu ~ dnorm(0, 0.000001)

# likelihood
for(i in 1:N){
p[i] <- 1/(sigma2[i])
y[i] ~ dnorm(theta[i], p[i])
theta[i] ~ dnorm (mu, tau)
y.rep[i] <- replicate.post(y[i])
}

## Posterior Checks
y.min <- ranked( y[1:N], 1 ) 
y.max <- ranked( y[1:N], N ) 
y.rep.min <- ranked( y.rep[1:N], 1 ) 
y.rep.max <- ranked( y.rep[1:N], N ) 
y.rep.min.dif <- y.min - y.rep.min 
y.rep.max.dif <- y.rep.max - y.max 

pred.y.min <- step( y.rep.min.dif ) 
pred.y.max <- step( y.rep.max.dif ) 

# Examine range 
y.rep.range.dif <- ( y.rep.max - y.rep.min ) - ( y.max - y.min ) 
pred.y.range <- step( y.rep.range.dif ) 

# Examine individual observations 
for( i in 1:N ) { 
y.rep.dif[i] <- y.rep[i] - y[i] 
pred.y[i] <- step( y.rep.dif[i] )
} 

# Examine sample and replicate sample mean and standard deviation 

pred.mean <- step( mean( y.rep[1:N] ) - mean( y[1:N] ) ) 
pred.sd <- step( sd(y.rep[1:N] ) - sd( y[1:N] ) )

}"
cat(re_model_3,file="./modelFiles/efm_re_3.txt")
re_sim_3 = bugs(
  data = list(N = nrow(efm_data[efm_data$study_type==1,]), y = efm_data$RD[efm_data$study_type==1],  sigma2 = efm_data$var_RD[efm_data$study_type==1]),
  inits = list(list(sd = sqrt(1/0.2), mu = -10),
               list(sd = sqrt(1/10), mu = 0),
               list(sd =  sqrt(1/20), mu = 10)),
  model.file = "C:/Users/AndyPohl/Dropbox/PhD - Calgary/Course_work/STAT619/Project/Analysis/EFM/modelFiles/efm_re_3.txt",
  parameters.to.save = c("tau", "var", "sd", "mu", "theta", "pred.mean", "pred.sd", "pred.y", "pred.y.max", "pred.y.min", "pred.y.range"),
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = burnin,
  n.thin =1,
  bugs.directory = 'c:/Program Files (x86)/OpenBUGS/OpenBUGS323/',
  working.directory = getwd(),
  program = 'OpenBUGS')
save(re_sim_3, file = './MCMC/re_sim_3.RData')


##### Mdl 4: Half Cauchy #####
re_model_4 = "
model{
# prior
sd ~ dt(0, 25, 1)I(0,)
var <- pow(sd, 2)
tau <- 1/var

mu ~ dnorm(0, 0.000001)

# likelihood
for(i in 1:N){
p[i] <- 1/(sigma2[i])
y[i] ~ dnorm(theta[i], p[i])
theta[i] ~ dnorm (mu, tau)
y.rep[i] <- replicate.post(y[i])
}

## Posterior Checks
y.min <- ranked( y[1:N], 1 ) 
y.max <- ranked( y[1:N], N ) 
y.rep.min <- ranked( y.rep[1:N], 1 ) 
y.rep.max <- ranked( y.rep[1:N], N ) 
y.rep.min.dif <- y.min - y.rep.min 
y.rep.max.dif <- y.rep.max - y.max 

pred.y.min <- step( y.rep.min.dif ) 
pred.y.max <- step( y.rep.max.dif ) 

# Examine range 
y.rep.range.dif <- ( y.rep.max - y.rep.min ) - ( y.max - y.min ) 
pred.y.range <- step( y.rep.range.dif ) 

# Examine individual observations 
for( i in 1:N ) { 
y.rep.dif[i] <- y.rep[i] - y[i] 
pred.y[i] <- step( y.rep.dif[i] )
} 

# Examine sample and replicate sample mean and standard deviation 
pred.mean <- step( mean( y.rep[1:N] ) - mean( y[1:N] ) ) 
pred.sd <- step( sd(y.rep[1:N] ) - sd( y[1:N] ) )

}"
cat(re_model_4,file="./modelFiles/efm_re_4.txt")
re_sim_4 = bugs(
  data = list(N = nrow(efm_data[efm_data$study_type==1,]), y = efm_data$RD[efm_data$study_type==1],  sigma2 = efm_data$var_RD[efm_data$study_type==1]),
  inits = list(list(sd = sqrt(1/0.2), mu = -10),
               list(sd = sqrt(1/10), mu = 0),
               list(sd =  sqrt(1/20), mu = 10)),
  model.file = "C:/Users/AndyPohl/Dropbox/PhD - Calgary/Course_work/STAT619/Project/Analysis/EFM/modelFiles/efm_re_4.txt",
  parameters.to.save = c("tau", "var", "sd", "mu", "theta", "pred.mean", "pred.sd", "pred.y", "pred.y.max", "pred.y.min", "pred.y.range"),
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = burnin,
  n.thin =1,
  bugs.directory = 'c:/Program Files (x86)/OpenBUGS/OpenBUGS323/',
  working.directory = getwd(),
  program = 'OpenBUGS')
save(re_sim_4, file = './MCMC/re_sim_4.RData')


##### Mdl 5: Uniform #####
re_model_5 = "
model{
# prior
sd ~ dunif(0,1000)
var <- pow(sd, 2)
tau <- 1/var

mu ~ dnorm(0, 0.000001)

# likelihood
for(i in 1:N){
p[i] <- 1/(sigma2[i])
y[i] ~ dnorm(theta[i], p[i])
theta[i] ~ dnorm (mu, tau)
y.rep[i] <- replicate.post(y[i])
}

## Posterior Checks
y.min <- ranked( y[1:N], 1 ) 
y.max <- ranked( y[1:N], N ) 
y.rep.min <- ranked( y.rep[1:N], 1 ) 
y.rep.max <- ranked( y.rep[1:N], N ) 
y.rep.min.dif <- y.min - y.rep.min 
y.rep.max.dif <- y.rep.max - y.max 

pred.y.min <- step( y.rep.min.dif ) 
pred.y.max <- step( y.rep.max.dif ) 

# Examine range 
y.rep.range.dif <- ( y.rep.max - y.rep.min ) - ( y.max - y.min ) 
pred.y.range <- step( y.rep.range.dif ) 

# Examine individual observations 
for( i in 1:N ) { 
y.rep.dif[i] <- y.rep[i] - y[i] 
pred.y[i] <- step( y.rep.dif[i] )
} 

# Examine sample and replicate sample mean and standard deviation 
pred.mean <- step( mean( y.rep[1:N] ) - mean( y[1:N] ) ) 
pred.sd <- step( sd(y.rep[1:N] ) - sd( y[1:N] ) )

}"
cat(re_model_5,file="./modelFiles/efm_re_5.txt")
re_sim_5 = bugs(
  data = list(N = nrow(efm_data[efm_data$study_type==1,]), y = efm_data$RD[efm_data$study_type==1],  sigma2 = efm_data$var_RD[efm_data$study_type==1]),
  inits = list(list(sd = sqrt(1/0.2), mu = -10),
               list(sd = sqrt(1/10), mu = 0),
               list(sd =  sqrt(1/20), mu = 10)),
  model.file = "C:/Users/AndyPohl/Dropbox/PhD - Calgary/Course_work/STAT619/Project/Analysis/EFM/modelFiles/efm_re_5.txt",
  parameters.to.save = c("tau", "var", "sd", "mu", "theta", "pred.mean", "pred.sd", "pred.y", "pred.y.max", "pred.y.min", "pred.y.range"),
  n.chains = n.chains,
  n.iter = n.iter,
  n.burnin = burnin,
  n.thin =1,
  bugs.directory = 'c:/Program Files (x86)/OpenBUGS/OpenBUGS323/',
  working.directory = getwd(),
  program = 'OpenBUGS')
save(re_sim_5, file = './MCMC/re_sim_5.RData')

### Including Observational Evidence

#independent analysis...
obs_model = "
model{
# prior
tau ~ dgamma(0.001, 0.001)
var <- 1/tau
sd <- sqrt(var)

mu ~ dnorm(0, 0.000001)


# likelihood
for(i in 1:N){
p[i] <- 1/(sigma2[i])
y[i] ~ dnorm(theta[i], p[i])
theta[i] ~ dnorm (mu, tau)
y.rep[i] <- replicate.post(y[i])
}

# Posterior Checks
y.min <- ranked( y[1:N], 1 ) 
y.max <- ranked( y[1:N], N ) 
y.rep.min <- ranked( y.rep[1:N], 1 ) 
y.rep.max <- ranked( y.rep[1:N], N ) 
y.rep.min.dif <- y.min - y.rep.min 
y.rep.max.dif <- y.rep.max - y.max 

pred.y.min <- step( y.rep.min.dif ) 
pred.y.max <- step( y.rep.max.dif ) 

# Examine range 
y.rep.range.dif <- ( y.rep.max - y.rep.min ) - ( y.max - y.min ) 
pred.y.range <- step( y.rep.range.dif ) 

# Examine individual observations 
for( i in 1:N ) { 
y.rep.dif[i] <- y.rep[i] - y[i] 
pred.y[i] <- step( y.rep.dif[i] )
} 

# Examine sample and replicate sample mean and standard deviation 
pred.mean <- step( mean( y.rep[1:N] ) - mean( y[1:N] ) ) 
pred.sd <- step( sd(y.rep[1:N] ) - sd( y[1:N] ) )
}"

cat(obs_model,file="./modelFiles/efm_obs.txt")


