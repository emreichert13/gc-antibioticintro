
# clean output of ODE solver into long format
SI.clean <- function(data){
  
data_long <- data %>%
  gather(key = type, value = individuals, -time)

#Classify infections by sexual risk group
data_long <- data_long %>%
  mutate(RiskGroup = ifelse(type == "S1" | type == "Y01" | type == "Ya1" | type == "Yb1" | type == "Yab1" | type == "Z01" | type == "Za1" | type == "Zb1" | type == "Zab1", "Low",
                            ifelse(type == "S2" | type == "Y02" | type == "Ya2" | type == "Yb2" | type == "Yab2" | type == "Z02" | type == "Za2" | type == "Zb2" | type == "Zab2", "Intermediate", 
                                   ifelse(type == "S3" | type == "Y03" | type == "Ya3" | type == "Yb3" | type == "Yab3" | type == "Z03" | type == "Za3" | type == "Zb3" | type == "Zab3", "High", NA))))
#classify infections by resistance profile & symp/asymp status
data_long <- data_long %>%
  mutate(Profile = str_sub(type, end=-2),
         InfectState = ifelse(grepl("S", type, fixed = TRUE), "Susceptible",
                              ifelse(grepl("Y", type, fixed = TRUE), "Sympomatic Infected",
                                     ifelse(grepl("Z", type, fixed = TRUE), "Asymptomatic Infected", NA))),
         ResistState = ifelse(Profile == "Y0" | Profile == "Z0", 'Neither',
                              ifelse(Profile == "Ya" | Profile == "Za", "A only",
                                     ifelse(Profile == "Yb" | Profile == "Zb", "B only",
                                            ifelse(Profile == "Yab" | Profile == "Zab", "A and B", NA)))))
return(data_long)
}

# quick visualization of results of ODE solver by Resistance Profile, sexual activity group, Infection status
SI.visualize <- function(data) { 
  
p1 <- data %>%
  group_by(time, ResistState) %>%
  summarise(N = sum(individuals)) %>%
  ggplot() + 
  geom_line(aes(x = time, y = log10(N), col = ResistState), size = 1.5) + 
  theme_light() + 
  ggtitle("GC Infections by Resistance Profile over Time") + ylim(c(0,6))

p2 <- data %>%
  group_by(time, ResistState, RiskGroup) %>%
  summarise(N = sum(individuals)) %>%
  ggplot() + 
  geom_line(aes(x = time, y = log10(N), col = ResistState), size = 1.5) + 
  theme_light() + 
  ggtitle("GC Infections by Resistance Profile over Time, by Risk") + ylim(c(0,6)) +
  facet_wrap(~RiskGroup, nrow = 3)
  

p3 <- data %>% 
  group_by(time, InfectState) %>% 
  summarise(N = sum(individuals)) %>%
  ggplot() + 
  geom_line(aes(x = time, y = log10(N), col = InfectState), size = 1.5) + 
  theme_light() + 
  scale_color_manual(values = c("deeppink3","midnightblue", "orangered")) +
  ggtitle("GC Infections by Clinical Profile over Time")

p4 <- data %>%
  group_by(time, InfectState, RiskGroup) %>%
  summarise(N = sum(individuals)) %>%
  ggplot() +
  geom_line(aes(x = time, y = log10(N), col = InfectState), size = 1.5) + 
  theme_light() + 
  facet_wrap(~RiskGroup) + 
  scale_color_manual(values = c("deeppink3","midnightblue", "orangered")) +
  ggtitle("GC Infections by Clinical Profile over Time, by Risk")
list(p1, p2, p3, p4)
}

# calculate summary stats from ODE output
# including time to loss (5% resistance) of Drug A, drug B, dual resistance
# as well as prevalence of gonorrhea at equilibrium, and @ time of loss
SI.outputs <- function(data, data_long) {
  end <- data[nrow(data),]
  #overall prevalence of gonorrhea at t=end
  prev_GC <- sum(end[,5:28])/sum(end[,2:28])
  
  #proportion of cases w/ some resistance at t=end
  prev_resistance <- sum(end[,11:28])/(sum(end[,5:28]))

  data_cond <- data_long %>%
    group_by(time, Profile) %>%
    summarise(N = sum(individuals)) %>%
    spread(., key = Profile, value = N) %>%
    mutate(prevA = round((Ya + Za + Yab + Zab)/(Y0 + Z0 + Ya + Za + Yb + Zb + Yab + Zab),3),
           prevB = round((Yb + Zb + Yab + Zab)/(Y0 + Z0 + Ya + Za + Yb + Zb + Yab + Zab),3),
           prevAB = round((Yab + Zab)/(Y0 + Z0 + Ya + Za + Yb + Zb + Yab + Zab),3))
  
  #Time to loss of one or both drugs (>5% prevalence of resistance)
  LossA_5 = min(data_cond$time[data_cond$prevA >= 0.05])/365
  LossB_5 = min(data_cond$time[data_cond$prevB >= 0.05])/365
  LossAB_5 = min(data_cond$time[data_cond$prevAB >= 0.05])/365
  
  #Time to loss of one or both drugs (>1% prevalence of resistance)
  LossA_1 = min(data_cond$time[data_cond$prevA >= 0.01])/365
  LossB_1 = min(data_cond$time[data_cond$prevB >= 0.01])/365
  LossAB_1 = min(data_cond$time[data_cond$prevAB >= 0.01])/365
  
  Loss_Each5 = max(LossA_5[is.finite(LossA_5)], LossB_5[is.finite(LossB_5)])*365
  time_loss = data[data$time == Loss_Each5,]
  prev_GC_loss <- sum(time_loss[,5:28])/sum(time_loss[,2:28])
  
  list(PrevGC = prev_GC, PrevRes = prev_resistance, PrevGCLoss = prev_GC_loss,
       LossA_5 = LossA_5, LossB_5 = LossB_5, LossAB_5 = LossAB_5,
       LossA_1 = LossA_1, LossB_1 = LossB_1, LossAB_1 = LossAB_1)
}

# calculate GC prevalence based on current parameters, for calibration
pred_fun_er <- function(params){
  params2 <- list(pop=10^6, pop.p=c(0.3, 0.6, 0.1), rho = 1/(20*365), 
                  prA = 1, omega_a = 10^-8, fA = 0.98, 
                  pi_s = 0.90)
  parms <- c(params, params2)
  years = 2
  dt <- seq(0, 365*years, 1) 
  calibration_sim <- as.data.frame(ode(inits, dt, calibration.SI, parms = parms))   
  calibration_sim <- round(calibration_sim, 0)
  end <- calibration_sim[nrow(calibration_sim),]
  #Overall prevalence of gonorrhea at t=end years
  prev_GC <- sum(end[,5:16])/sum(end[,2:16])
  prev_GC
  return(prev_GC)
}

# estimate alpha and beta for our beta distribution of prevalence of GC, for calibration
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# Calibration Model (Drug A Only)
calibration.SI <- function(t, x, parms){
  with(as.list(c(t, x, parms)),{
    N = c(N1, N2, N3) #total population
    S = c(S1, S2, S3) #not infected
    Y0 = c(Y01, Y02, Y03) #symptomatic infected, no resistance
    Ya = c(Ya1, Ya2, Ya3) #symptomatic infected, resistant to A
    Z0 = c(Z01, Z02, Z03) #asymptomatic infected, no resistance
    Za = c(Za1, Za2, Za3) #asymptomatic infected, resistant to A
    #makes 3x3 contact matrix
    activities <- c(1*c_min/365, 5*c_min/365, 20*c_min/365)
    beta <- (1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + 
      epsilon*activities/(pop*pop.p)*diag(3)
    beta <- beta * b #contacts * transmission pr per partnership
    #susceptibles
    dS <- -beta %*% ((Y0 + Z0) + fA*(Ya+Za)) *S + 
      (1-omega_a*prA)*(Ts*Y0 + Tm*Z0) +
      pi_s*Tsr*Ya +
      g*(Y0 + Ya + Z0 + Za) + 
      rho*(S + Y0 + Ya + Z0 + Za) - rho*S
    #infections w/ no resistance
    dY0 <- sigma*(beta %*% (Y0+Z0)*S) - Ts*Y0 - g*Y0 - rho*Y0
    dZ0 <- (1-sigma)*(beta %*% (Y0+Z0)*S) - Tm*Z0 - g*Z0 - rho*Z0
    #infections w/ resistance to A
    dYa <- sigma*fA*(beta %*% (Ya+Za)*S) +
      omega_a*prA*Ts*Y0 - 
      pi_s*Tsr*Ya -
      g*Ya - rho*Ya
    dZa <- (1-sigma)*fA*(beta %*% (Ya+Za)*S) +
      omega_a*prA*Tm*Z0 -
      g*Za - rho*Za
    der <- c(dS, dY0, dZ0, dYa, dZa)
    list(der)
  })
}

logit<-function(x) {log(x/(1-x))}

ilogit <-function(x) {1/(1+exp(-x))}
