
#E Reichert 2022

#Run Sensitivity Analysis + Save Output
#Explore effect of parameter values for fB and omegaB

#UPDATE to your file path
knit("~/Documents/2021 Grad lab research/GCtransmission.Rmd")

# Set parameters that will stay constant throughout model
# NOTE - these parameters will produce results seen in Figs 2-3 (baseline params for drug A)
# Params for Drug A can be altered to produce output as seen in Supp. Fig. 4
pop = 10^6                        #pop size
pop.p = c(0.3, 0.6, 0.1)          #relative size of each risk group; lo, M, hi
c_min = 1.2165433
activities = c(1*c_min/365, 
               5*c_min/365, 
               20*c_min/365)      #sexual contacts per day
epsilon = 0.2427385               #mixing parameter

#PARAM VALUES
b = 0.4566313        #transmission Pr per partnership
sigma = 0.6009376    #Pr of symptomatic infection
g = 1/(168.5402978)  #natural recovery rate from infection
Ts =1/(11.4216522)   #time to treatment for symptomatic infection
Tm = 0.4035764/365   #screening rate (time to treatment for asymptomatic infection)
rho = 1/(20*365)     #model entry/exit rate
omega_a = 10^-8 #10^-4  #Pr of emergence of resistance on treatment with A (ceftriaxone)
#omega_b = 10^-4     #Pr of emergence of resistance on treatment with B (new drug)
fA = 0.98 #0.90      #relative fitness, resistant to A
#fB = 0.95           #relative fitness, resistant to B (is this reasonable??)
#fAB = fA*fB         #relative fitness, dual resistance
pi_s = 0.90          #Pr of retreatment if initial treatment failure, symptomatic
Tsr = Ts/3           #time to retreatment for symptomatic infection, if failure
resA = 0.0001        #initial prev of resistance to A

###############  LONG TERM results ########################

years = 100
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

#estimated proportion of symp. infections at cross-sectional point in time = 10%
prev_symp <- 0.10
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

# set initial conditions
inits <- c(S1 = N1-gc_lo, S2 = N2-gc_md, S3 = N3-gc_hi, 
           Y01 = Y0[1], Y02 = Y0[2], Y03 = Y0[3],
           Z01 = Z0[1], Z02 = Z0[2], Z03 = Z0[3], 
           Ya1 = YA[1], Ya2 = YA[2], Ya3 = YA[3], 
           Za1 = ZA[1], Za2 = ZA[2], Za3 = ZA[3], 
           Yb1 = 0, Yb2 = 0, Yb3 = 0, 
           Zb1 = 0, Zb2 = 0, Zb3 = 0, 
           Yab1 = 0, Yab2 = 0, Yab3 = 0, 
           Zab1 = 0, Zab2 = 0, Zab3 = 0) 

dt <- seq(0, 365*years, tstep)

#Set parameter values for fB and omegaB that we will explore
fB_range <- seq(0.80, 1, 0.05) #ranges from 0-20% fitness cost
omega_b_range <- c(10^seq(-10,-2)) #ranges for Pr of resistance emergence upon treatment

prA = 0.5            #Pr of treatment with A
prB = 0.5            #Pr of treatment with B

#run ODE models over all combinations of these params
SensAnalyze <- function (fB_range, omega_b_range) {
  #create data frame with 0 rows and 3 columns
  df <- data.frame(matrix(ncol = 7, nrow = 0))
  #provide column names
  colnames(df) <- c('fB', 'omega_B', 'PrevGC', 'LossA', 'LossB', 'LossAB', 'Strategy')
  for(i in fB_range)
  {
    for(j in omega_b_range)
    {
      fB <- i
      omega_b <- j
      fAB <- fA*fB
      #list out parameters defined above
      parms <- list(g=g, pop=pop, pop.p=pop.p, epsilon=epsilon, b = b, sigma = sigma, 
                    Ts = Ts, Tm = Tm, rho = rho, prA = prA, prB = prB, omega_a = omega_a, 
                    omega_b = omega_b, fA = fA, fB = fB, fAB = fAB, pi_s = pi_s, 
                    Tsr = Tsr)
      #Strategy 1
      parms$prA <- 0.5
      parms$prB <- 0.5
      random_sim <- as.data.frame(ode(inits, dt, random.SI, parms=parms)) 
      random_sim_long <- SI.clean(random_sim)
      outputs <- SI.outputs(random_sim, random_sim_long)
      results1 <- data.frame(fB = i, omega_B = j, PrevGC = outputs$PrevGC, LossA = outputs$LossA_5, 
                             LossB = outputs$LossB_5, LossAB = outputs$LossAB_5, Strategy = 1)
      #mono for calculation of Q
      parms$prA <- 1
      parms$prB <- 0
      mono_sim <- as.data.frame(ode(inits, dt, random.SI, parms=parms)) 
      mono_sim_long <- SI.clean(mono_sim)
      outputs <- SI.outputs(mono_sim, mono_sim_long)
      Q <- outputs$LossA_5*365 #time 5% resistance threshold met for drug A
      results0 <- data.frame(fB = i, omega_B = j, PrevGC = outputs$PrevGC, LossA = outputs$LossA_5, 
                             LossB = outputs$LossB_5, LossAB = outputs$LossAB_5, Strategy = 0)
      
      #Strategy 2
      combo_sim <- as.data.frame(ode(inits, dt, combo.SI, parms=parms)) 
      combo_sim_long <- SI.clean(combo_sim)
      outputs <- SI.outputs(combo_sim, combo_sim_long)
      results2 <- data.frame(fB = i, omega_B = j, PrevGC = outputs$PrevGC, LossA = outputs$LossA_5, 
                             LossB = outputs$LossB_5, LossAB = outputs$LossAB_5, Strategy = 2)
      #Strategy 3
      reserve_sim <- as.data.frame(ode(inits, dt, reserve.SI, parms=c(parms, Q = Q))) 
      reserve_sim_long <- SI.clean(reserve_sim)
      outputs <- SI.outputs(reserve_sim, reserve_sim_long)
      results3 <- data.frame(fB = i, omega_B = j, PrevGC = outputs$PrevGC, LossA = outputs$LossA_5, 
                             LossB = outputs$LossB_5, LossAB = outputs$LossAB_5, Strategy = 3)
      #Strategy 4
      gradual_sim <- as.data.frame(ode(inits, dt, gradual.SI, parms=parms)) 
      gradual_sim_long <- SI.clean(gradual_sim)
      outputs <- SI.outputs(gradual_sim, gradual_sim_long)
      results4 <- data.frame(fB = i, omega_B = j, PrevGC = outputs$PrevGC, LossA = outputs$LossA_5, 
                             LossB = outputs$LossB_5, LossAB = outputs$LossAB_5, Strategy = 4)
      
      df <- rbind(df, results0, results1, results2, results3, results4)
    }
  }
  return(df)
}

#save data to desired location
SensResults <- SensAnalyze(fB_range, omega_b_range)
write.csv(SensResults, "SensResults_rerun.csv")


###############  SHORT TERM results ########################

dt <- seq(0, 365*20, tstep)

SensAnalyze_short <- function (fB_range, omega_b_range) {
  #create data frame with 0 rows and 3 columns
  df <- data.frame(matrix(ncol = 7))
  colnames(df) <- colnames(random_sim_long)
  df <- df %>%
    mutate(fB = NA, omega_b = NA, Strategy = NA)
  #provide column names
  for(i in fB_range)
  {
    for(j in omega_b_range)
    {
      fB <- i
      omega_b <- j
      fAB <- fA*fB
      #list out parameters defined above
      parms <- list(g=g, pop=pop, pop.p=pop.p, epsilon=epsilon, b = b, sigma = sigma, 
                    Ts = Ts, Tm = Tm, rho = rho, prA = prA, prB = prB, omega_a = omega_a, 
                    omega_b = omega_b, fA = fA, fB = fB, fAB = fAB, pi_s = pi_s, 
                    Tsr = Tsr)
      #Strategy 1
      parms$prA <- 0.5
      parms$prB <- 0.5
      random_sim <- as.data.frame(ode(inits, dt, random.SI, parms=parms)) %>% 
        select(-prevA) %>%
        mutate(CumCases = cumsum(Inc))
      random_sim_long <- SI.clean(random_sim) %>%
        mutate(fB = fB, omega_b = omega_b, Strategy = 1) %>%
        filter(time == 3650 | time == 7300)
      #mono for calculation of Q
      parms$prA <- 1
      parms$prB <- 0
      mono_sim <- as.data.frame(ode(inits, dt, random.SI, parms=parms)) 
      mono_sim_long <- SI.clean(mono_sim)
      outputs <- SI.outputs(mono_sim, mono_sim_long)
      Q <- outputs$LossA_5*365
      #Strategy 2
      combo_sim <- as.data.frame(ode(inits, dt, combo.SI, parms=parms)) %>%
        mutate(CumCases = cumsum(Inc))
      combo_sim_long <- SI.clean(combo_sim) %>%
        mutate(fB = fB, omega_b = omega_b, Strategy = 2) %>%
        filter(time == 3650 | time == 7300)
      #Strategy 3
      reserve_sim <- as.data.frame(ode(inits, dt, reserve.SI, parms=c(parms, Q = Q))) %>%
        select(-c(prevA, prA, prB)) %>%
        mutate(CumCases = cumsum(Inc))
      reserve_sim_long <- SI.clean(reserve_sim) %>%
        mutate(fB = fB, omega_b = omega_b, Strategy = 3) %>%
        filter(time == 3650 | time == 7300)
      #Strategy 4
      gradual_sim <- as.data.frame(ode(inits, dt, gradual.SI, parms=parms)) %>%
        mutate(CumCases = cumsum(Inc))
      gradual_sim_long <- SI.clean(gradual_sim) %>%
        mutate(fB = fB, omega_b = omega_b, Strategy = 4) %>%
        filter(time == 3650 | time == 7300)
      
      df <- rbind(df, random_sim_long, combo_sim_long, reserve_sim_long, gradual_sim_long)
    }
  }
  return(df)
}

#save data to desired location
SensResults_new <- SensAnalyze_short(fB_range, omega_b_range)
write.csv(SensResults_new, "SensResults_new_update.csv")

