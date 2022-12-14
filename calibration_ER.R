
#Modeling various antibiotic introduction strategies for N. gonorrhoeae
#E Reichert, 04.2022

#Model Calibration
#see what parameters lead to under status quo conditions

source("~/Documents/2021 Grad lab research/GCfunctions.R")
library(deSolve)

#SET PARAM VALUES
#parameters that will stay constant throughout model
pop = 10^6                                           #pop size
pop.p = c(0.3, 0.6, 0.1)                             #relative size of each risk group; low, M, high
#activities = c(1*1.16/365, 5*1.16/365, 20*1.16/365)  #sexual contacts per day

rho = 1/(20*365)     #model entry/exit rate
omega_a = 10^-8      #Pr of emergence of resistance on treatment with A (ceftriaxone)
#omega_b = 0.01      #Pr of emergence of resistance on treatment with B (new drug)
prA = 1              #Pr of treatment with A
#prB = 0.5           #Pr of treatment with B
fA = 0.98            #relative fitness, resistant to A
#fB = 0.90           #relative fitness, resistant to B (is this reasonable??)
#fAB = fA*fB         #relative fitness, dual resistance
pi_s = 0.9          #Pr of retreatment if initial treatment failure, symptomatic
resA = 0.0001        #initial prev of resistance to A

#Set model duration + initial conditions

#how long to run the model?
years = 2
tstep = 1 #in days

#set N for each sexual risk group (1 = low risk, 2 = med risk, 3 = hi risk)
N <- pop
N1 <- N*pop.p[1]
N2 <- N*pop.p[2]
N3 <- N*pop.p[3]

#distribute GC cases to have overall 3% prevalence
x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)

gc_lo <- round(N1*x*0.029,0)
gc_md <- round(N2*x*0.154,0)
gc_hi <- round(N3*x*0.817,0)

#estimated proportion of symp. infections at cross-sectional point in time
prev_symp <- 0.089
Y_gc_lo <- gc_lo*prev_symp
Z_gc_lo <- gc_lo*(1-prev_symp)
Y_gc_md <- gc_md*prev_symp
Z_gc_md <- gc_md*(1-prev_symp)
Y_gc_hi <- gc_hi*prev_symp
Z_gc_hi <- gc_hi*(1-prev_symp)

Y0 <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * (1-resA)
Z0 <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * (1-resA)

YA <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resA
ZA <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resA

inits <- c(S1 = N1-gc_lo, S2 = N2-gc_md, S3 = N3-gc_hi, 
           Y01 = Y0[1], Y02 = Y0[2], Y03 = Y0[3],
           Z01 = Z0[1], Z02 = Z0[2], Z03 = Z0[3], 
           Ya1 = YA[1], Ya2 = YA[2], Ya3 = YA[3], 
           Za1 = ZA[1], Za2 = ZA[2], Za3 = ZA[3]) 

#create vector w/ all timepoints
dt <- seq(0, 365*years, tstep) 

##########################################
### model parameters from MLE fitting  ###
##########################################

#start with parameter values from Tuite et al. 2017
theta <- c(logit.b=logit(0.44),  #transmission probability
           log.c_min=log(1.16),  #min rate of partner change
           logit.epsilon=logit(0.23), #epsilon
           logit.sigma = logit(0.60), #pr of incident symptomatic infection
           log.Ts=log(0.026*365), #duration of infectiousness if symptomatic (days) = average time to treatment
           log.g=log(.5*365), #duration of infectiousness if asymptomatic and untreated (days)
           logit.Tm =logit(0.39)) #screening rate per year


### set parameters ###
b <- ilogit(theta["logit.b"])
c_min <- exp(theta["log.c_min"]) 
epsilon <- ilogit(theta["logit.epsilon"])
sigma <- ilogit(theta["logit.sigma"])
Ts <- exp(theta["log.Ts"])/365
g <- exp(theta["log.g"])/365
Tm <- ilogit(theta["logit.Tm"])


#######################################
#### Maximum likelihood estimation ####
#######################################
#used for estimating parameters

model.epi.loglik <- function(theta) {
  b <- ilogit(theta["logit.b"])
  c_min <- exp(theta["log.c_min"]) 
  epsilon <- ilogit(theta["logit.epsilon"])
  sigma <- ilogit(theta["logit.sigma"])
  Ts <- 1/(exp(theta["log.Ts"]))
  g <- 1/(exp(theta["log.g"]))
  Tm <- ilogit(theta["logit.Tm"])/365
  params <-list(c_min = c_min, epsilon = epsilon, sigma = sigma, b=b,Ts = Ts,Tm = Tm, g = g, Tsr = Ts/3)
  pred <- pred_fun_er(params)
  beta.params.prev <- estBetaParams(mu = pred, var = 1.47e-5)
  ll <- sum(dbeta(x= 0.03, beta.params.prev$alpha, beta.params.prev$beta, log=TRUE)) #calculate likelihood
  ll[is.na(ll)]<-(-1e20)
  print(ll)
  c(prev=pred,ll=ll)
}


f.optim <-function(theta) {
  Res <-  model.epi.loglik(theta)
  LogLL <- Res["ll"] #Model is returning the Log likelihood
  return(-LogLL)
}

#### run MLE optimization ###
library(bbmle)
values.start=theta
parnames(f.optim)<-names(values.start)
fit0 <- bbmle::mle2(f.optim, start=values.start,  vecpar=TRUE,  optimizer="optim"); fit0
fit <- mle2(f.optim, start=coef(fit0),  vecpar=TRUE, optimizer="optim"); fit
theta.fit<-coef(fit)

exp(theta.fit)
ilogit(theta.fit)
################################################

#Check calibration model fit with these params

#SET PARAM VALUES
#parameters that will stay constant throughout model
pop = 10^6                      #pop size
pop.p = c(0.3, 0.6, 0.1)        #relative size of each risk group; low, M, high
c_min = 1.2165433
#activities = c(1*c_min/365, 5*c_min/365, 20*c_min/365)  #sexual contacts per day
epsilon = 0.2427385                                      #mixing parameter

b = 0.4566313        #transmission Pr per partnership
sigma = 0.6009376    #Pr of symptomatic infection
g = 1/(168.5402978)  #natural recovery rate from infection
Ts =1/(11.4216522)   #time to treatment for symptomatic infection
Tm = 0.4035764/365   #screening rate (time to treatment for asymptomatic infection)
rho = 1/(20*365)     #model entry/exit rate
omega_a = 10^-8      #Pr of emergence of resistance on treatment with A (ceftriaxone)
#omega_b = 0.01      #Pr of emergence of resistance on treatment with B (new drug)
prA = 1              #Pr of treatment with A
#prB = 0.5           #Pr of treatment with B
fA = 0.98            #relative fitness, resistant to A
#fB = 0.90           #relative fitness, resistant to B (is this reasonable??)
#fAB = fA*fB         #relative fitness, dual resistance
pi_s = 0.90          #Pr of retreatment if initial treatment failure, symptomatic
Tsr = Ts/3           #time to retreatment for symptomatic infection, if failure
resA = 0.0001        #initial prev of resistance to A

#Set model duration + initial conditions
years = 2
tstep = 1 #in days

#set N for each sexual risk group (1 = low risk, 2 = med risk, 3 = hi risk)
N <- pop
N1 <- N*pop.p[1]
N2 <- N*pop.p[2]
N3 <- N*pop.p[3]

#distribute GC cases to have overall 3% prevalence
x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)

gc_lo <- round(N1*x*0.029,0)
gc_md <- round(N2*x*0.154,0)
gc_hi <- round(N3*x*0.817,0)

#estimated proportion of symp. infections at cross-sectional point in time
prev_symp <- 0.10 #(.6*(.026))/(.6*(.026)+.4*(.39))
Y_gc_lo <- gc_lo*prev_symp
Z_gc_lo <- gc_lo*(1-prev_symp)
Y_gc_md <- gc_md*prev_symp
Z_gc_md <- gc_md*(1-prev_symp)
Y_gc_hi <- gc_hi*prev_symp
Z_gc_hi <- gc_hi*(1-prev_symp)

Y0 <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * (1-resA)
Z0 <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * (1-resA)

YA <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resA
ZA <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resA

#start with overall prevalence of 3% gonorrhea
inits <- c(S1 = N1-gc_lo, S2 = N2-gc_md, S3 = N3-gc_hi, 
           Y01 = Y0[1], Y02 = Y0[2], Y03 = Y0[3],
           Z01 = Z0[1], Z02 = Z0[2], Z03 = Z0[3], 
           Ya1 = YA[1], Ya2 = YA[2], Ya3 = YA[3], 
           Za1 = ZA[1], Za2 = ZA[2], Za3 = ZA[3]) 

#create vector w/ all timepoints
dt <- seq(0, 365*years, tstep) 

#list out parameters
parms <- list(g=g, pop=pop, pop.p=pop.p, epsilon=epsilon, b = b, 
              sigma = sigma, Ts = Ts, Tm = Tm, rho = rho, prA = prA, 
              omega_a = omega_a, fA = fA, pi_s = pi_s, Tsr = Tsr, c_min = c_min)

#Run the model
calibration_sim <- as.data.frame(ode(inits, dt, calibration.SI, parms=parms))   

#ensure that total pop size does not change
total_check <- calibration_sim %>% mutate(RowSum=rowSums(.[setdiff(names(.),"time")])) %>% select(RowSum)

#Data tidying
calibration_sim_long <- SI.clean(calibration_sim)

# Visualize some outputs
SI.visualize(data = calibration_sim_long)

end <- calibration_sim[nrow(calibration_sim),]

#Overall prevalence of gonorrhea at t=end years
prev_GC_calibration <- sum(end[,5:16])/sum(end[,2:16])
prev_GC_calibration

#overall prevalence of resistance among cases at t=end years
prev_resistance_calibration <- sum(end[,11:16])/(sum(end[,5:16]))
prev_resistance_calibration

######## APPENDIX -- explore params for prevalence distribution
#define range
p = seq(0, 0.10, length=100)
#create plot of Beta distribution with shape parameters 2 and 10
vars <- estBetaParams(mu = 0.03, var = 1.47e-5)
plot(p, dbeta(p, vars$alpha, vars$beta), type='l')
plot(p, pbeta(p, vars$alpha, vars$beta), type='l')

# Sample size
n = 1000
# Parameters of the beta distribution
alpha = vars$alpha
beta = vars$beta
# Simulate some data
set.seed(1)
x = rbeta(n, alpha, beta)

# Note that the distribution is not symmetrical
curve(dbeta(x,alpha,beta))

vars$alpha/(vars$alpha + vars$beta)

# Negative log likelihood for the beta distribution
nloglikbeta = function(mu, sig) {
  alpha = mu^2*(1-mu)/sig^2-mu
  beta = alpha*(1/mu-1)
  -sum(dbeta(x, alpha, beta, log=TRUE))
}

library(stats4)
est = mle(nloglikbeta, start=list(mu=mean(x), sig=sd(x)))
confint(est)
t.test(x)$conf.int

