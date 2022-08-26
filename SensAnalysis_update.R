#Sensitivity Analysis
#Explore effect of values of fB and omegaB

library(viridis)
library(gridExtra)

source("~/Documents/2021 Grad lab research/GCtransmission.Rmd")

# Set parameters that will stay constant throughout model
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
omega_a = 10^-8      #Pr of emergence of resistance on treatment with A (ceftriaxone)
#omega_b = 10^-4     #Pr of emergence of resistance on treatment with B (new drug)
fA = 0.98            #relative fitness, resistant to A
#fB = 0.95           #relative fitness, resistant to B (is this reasonable??)
#fAB = fA*fB         #relative fitness, dual resistance
pi_s = 0.90          #Pr of retreatment if initial treatment failure, symptomatic
Tsr = Ts/3           #time to retreatment for symptomatic infection, if failure
resA = 0.0001        #initial prev of resistance to A

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

SensResults <- SensAnalyze(fB_range, omega_b_range)

#SensResults <- SensResults_fA
write.csv(SensResults, "SensResults_omegaA.csv")


################ ANALYZE RESULTS ####################

sens_baseline <- read.csv("SensResults_Baseline.csv") %>% select(-X)
sens_fA <- read.csv("SensResults_fA.csv") %>% select(-X)
sens_omegaA <- read.csv("SensResults_omegaA.csv") %>% select(-X)

SensResults <- sens_baseline

SensResults$LossA[is.infinite(SensResults$LossA) | SensResults$LossA > 100] <- 100
SensResults$LossB[is.infinite(SensResults$LossB)| SensResults$LossB > 100] <- 100
SensResults$LossAB[is.infinite(SensResults$LossAB) | SensResults$LossAB > 100] <- 100
SensResults <- SensResults %>%
  filter(Strategy != 0)

SensResults$LossBoth <- apply(SensResults[, 4:5], 1, max)
SensResults$LossOne <- apply(SensResults[, 4:5], 1, min)

SensResults$Strategy[SensResults$Strategy == 1] <- "1. 50-50"
SensResults$Strategy[SensResults$Strategy == 2] <- "2. Combination"
SensResults$Strategy[SensResults$Strategy == 3] <- "3. Reserve"
SensResults$Strategy[SensResults$Strategy == 4] <- "4. Gradual"

check <- SensResults %>%
  dplyr::group_by(fB, omega_B) %>%
  summarise(Strategy = Strategy[which.max(LossBoth)], value = max(LossBoth)) %>%
  mutate(Optimal = 1) %>%
  filter(value != 100) 

# SensResults$LossBoth[SensResults$LossBoth == 100] <- NA
# SensResults$LossOne[SensResults$LossOne == 100] <- NA

# calculate results relative to the reserve strategy
reserve <- SensResults %>% filter(Strategy == "3. Reserve") %>%
  select(fB, omega_B, LossBoth_reserve = LossBoth, LossOne_reserve = LossOne)
 
SensResults <- left_join(SensResults, reserve, by = c("fB", "omega_B"))
 
SensResults <- SensResults %>%
  mutate(Diff_both = LossBoth - LossBoth_reserve,
         Diff_one = LossOne - LossOne_reserve,
         Fold_both = LossBoth/LossBoth_reserve) %>% 
  filter(Strategy != "3. Reserve")

SensResults <- SensResults %>%
  filter(!(fB == 0.80 & omega_B <= 10^-6))

SensResults$fB[SensResults$fB == 1] <- "fB = 1"
SensResults$fB[SensResults$fB == 0.95] <- "fB = 0.95"
SensResults$fB[SensResults$fB == 0.90] <- "fB = 0.90"
SensResults$fB[SensResults$fB == 0.85] <- "fB = 0.85"
SensResults$fB[SensResults$fB == 0.80] <- "fB = 0.80"

SensResults <- SensResults %>%
  mutate(Label = ifelse(LossBoth == 100, 
                        paste(">", round(Diff_both,0), sep = ""),
                        round(Diff_both,0)),
         Diff_both = ifelse(Diff_both > 60, 60, Diff_both))

library(paletteer)
pdf("Figure2.pdf", width = 7.5, height = 3.5)
#jpeg("Figure2.jpg", width = 7.5, height = 3.5, units = "in", res = 300)
ggplot(SensResults, aes(Strategy, reorder(omega_B, desc(omega_B)), fill= Diff_both)) + 
  geom_tile() + theme_classic() +
  geom_text(aes(label = Label), col = "white", size = 2.5) +
  xlab("Strategy") + 
  ylab("Probability of de novo Resistance, Drug B") + 
  labs(fill = "Additional Time to Loss \nof Both Drugs, relative to \nReserve Strategy (years)") + 
  scale_fill_paletteer_c("grDevices::SunsetDark", direction = -1, na.value = "white") +
  theme(text = element_text(size=9), axis.text.x = element_text(angle = -45, hjust = 0.1, vjust = 0.1),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=7),
        strip.background = element_rect(fill = "deeppink3", colour = "white"),
        strip.text = element_text(color = "white", size = 8),
        axis.text=element_text(size=7),
        axis.title=element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 8)) +
  ggtitle("fB = Relative Fitness of Drug B Resistant Strains") +
  facet_grid(~fB)
dev.off()

#now, for altered drug A fitness
sens_fA$LossA[is.infinite(sens_fA$LossA) | sens_fA$LossA > 100] <- 100
sens_fA$LossB[is.infinite(sens_fA$LossB)| sens_fA$LossB > 100] <- 100
sens_fA$LossAB[is.infinite(sens_fA$LossAB) | sens_fA$LossAB > 100] <- 100
sens_fA <- sens_fA %>%
  filter(Strategy != 0 & Strategy != 4.2)

sens_fA$LossBoth <- apply(sens_fA[, 4:5], 1, max)
sens_fA$LossOne <- apply(sens_fA[, 4:5], 1, min)

sens_fA$Strategy[sens_fA$Strategy == 1] <- "1. 50-50"
sens_fA$Strategy[sens_fA$Strategy == 2] <- "2. Combination"
sens_fA$Strategy[sens_fA$Strategy == 3] <- "3. Reserve"
sens_fA$Strategy[sens_fA$Strategy == 4.1] <- "4. Gradual"

check <- sens_fA %>%
  dplyr::group_by(fB, omega_B) %>%
  summarise(Strategy = Strategy[which.max(LossBoth)], value = max(LossBoth)) %>%
  mutate(Optimal = 1) %>%
  filter(value != 100) 

# sens_fA$LossBoth[sens_fA$LossBoth == 100] <- NA
# sens_fA$LossOne[sens_fA$LossOne == 100] <- NA

# calculate results relative to the Reserve Strategy
reserve <- sens_fA %>% filter(Strategy == "3. Reserve") %>%
  select(fB, omega_B, LossBoth_reserve = LossBoth, LossOne_reserve = LossOne)

sens_fA <- left_join(sens_fA, reserve, by = c("fB", "omega_B"))
 
sens_fA <- sens_fA %>%
  mutate(Diff_both = LossBoth - LossBoth_reserve,
         Diff_one = LossOne - LossOne_reserve,
         Fold_both = LossBoth/LossBoth_reserve) %>% 
  filter(Strategy != "3. Reserve")

sens_fA <- sens_fA %>%
  mutate(Diff_both = ifelse(fB == 0.80 | (fB == 0.85 & omega_B <= 10^-4), NA, Diff_both))

sens_fA$fB[sens_fA$fB == 1] <- "fB = 1"
sens_fA$fB[sens_fA$fB == 0.95] <- "fB = 0.95"
sens_fA$fB[sens_fA$fB == 0.90] <- "fB = 0.90"
sens_fA$fB[sens_fA$fB == 0.85] <- "fB = 0.85"
sens_fA$fB[sens_fA$fB == 0.80] <- "fB = 0.80"

sens_fA <- sens_fA %>%
  mutate(Label = ifelse(LossBoth == 100, 
                        paste(">", round(Diff_both,0), sep = ""),
                        round(Diff_both,0)),
         Diff_both = ifelse(Diff_both > 80, 80, Diff_both))

g1 <- ggplot(sens_fA, aes(Strategy, reorder(omega_B, desc(omega_B)), fill= Diff_both)) + 
  geom_tile() + theme_classic() +
  geom_text(aes(label = Label), col = "white", size = 2.5) +
  xlab("Strategy") + 
  ylab("Probability of de novo Resistance, Drug B") + 
  labs(fill = "Additional Time to Loss \nof Both Drugs, relative to \nReserve Strategy (years)",
       title = "A. Relative Fitness of Drug A Resistant Strains (fA) = 0.90",
       subtitle = "fB = Relative Fitness of Drug B Resistant Strains") + 
  scale_fill_paletteer_c("grDevices::SunsetDark", direction = -1, na.value = "white") +
  theme(text = element_text(size=9), axis.text.x = element_text(angle = -45, hjust = 0.1, vjust = 0.1),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=7),
        strip.background = element_rect(fill = "deeppink3", colour = "white"),
        strip.text = element_text(color = "white", size = 8),
        axis.text=element_text(size=7),
        axis.title=element_text(size = 8),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        plot.title = element_text(size = 10)) +
  facet_grid(~fB)

#now, for altered omega A
sens_omegaA$LossA[is.infinite(sens_omegaA$LossA) | sens_omegaA$LossA > 100] <- 100
sens_omegaA$LossB[is.infinite(sens_omegaA$LossB)| sens_omegaA$LossB > 100] <- 100
sens_omegaA$LossAB[is.infinite(sens_omegaA$LossAB) | sens_omegaA$LossAB > 100] <- 100
sens_omegaA <- sens_omegaA %>%
  filter(Strategy != 0 & Strategy != 4.2)

sens_omegaA$LossBoth <- apply(sens_omegaA[, 4:5], 1, max)
sens_omegaA$LossOne <- apply(sens_omegaA[, 4:5], 1, min)

sens_omegaA$Strategy[sens_omegaA$Strategy == 1] <- "1. 50-50"
sens_omegaA$Strategy[sens_omegaA$Strategy == 2] <- "2. Combination"
sens_omegaA$Strategy[sens_omegaA$Strategy == 3] <- "3. Reserve"
sens_omegaA$Strategy[sens_omegaA$Strategy == 4.1] <- "4. Gradual"

check <- sens_omegaA %>%
  dplyr::group_by(fB, omega_B) %>%
  summarise(Strategy = Strategy[which.max(LossBoth)], value = max(LossBoth)) %>%
  mutate(Optimal = 1) %>%
  filter(value != 100) 

# sens_omegaA$LossBoth[sens_omegaA$LossBoth == 100] <- NA
# sens_omegaA$LossOne[sens_omegaA$LossOne == 100] <- NA

# calculate results relative to the Reserve Strategy
reserve <- sens_omegaA %>% filter(Strategy == "3. Reserve") %>%
  select(fB, omega_B, LossBoth_reserve = LossBoth, LossOne_reserve = LossOne)
 
sens_omegaA <- left_join(sens_omegaA, reserve, by = c("fB", "omega_B"))
 
sens_omegaA <- sens_omegaA %>%
  mutate(Diff_both = LossBoth - LossBoth_reserve,
         Diff_one = LossOne - LossOne_reserve,
         Fold_both = LossBoth/LossBoth_reserve) %>% 
  filter(Strategy != "3. Reserve")

sens_omegaA <- sens_omegaA %>%
  mutate(Diff_both = ifelse(fB == 0.80 & omega_B <= 10^-7, NA, Diff_both))

sens_omegaA$fB[sens_omegaA$fB == 1] <- "fB = 1"
sens_omegaA$fB[sens_omegaA$fB == 0.95] <- "fB = 0.95"
sens_omegaA$fB[sens_omegaA$fB == 0.90] <- "fB = 0.90"
sens_omegaA$fB[sens_omegaA$fB == 0.85] <- "fB = 0.85"
sens_omegaA$fB[sens_omegaA$fB == 0.80] <- "fB = 0.80"

sens_omegaA <- sens_omegaA %>%
  mutate(Label = ifelse(LossBoth == 100, 
                        paste(">", round(Diff_both,0), sep = ""),
                        round(Diff_both,0)),
         Diff_both = ifelse(Diff_both > 60, 60, Diff_both))

g2 <- ggplot(sens_omegaA, aes(Strategy, reorder(omega_B, desc(omega_B)), fill= Diff_both)) + 
  geom_tile() + theme_classic() +
  geom_text(aes(label = Label), col = "white", size = 2.5) +
  xlab("Strategy") + 
  ylab("Probability of de novo Resistance, Drug B") + 
  labs(fill = "Additional Time to Loss \nof Both Drugs, relative to \nReserve Strategy (years)",
       title = "B. Probability of de novo Resistance, Drug A = 1e-04",
       subtitle = "fB = Relative Fitness of Drug B Resistant Strains") + 
  scale_fill_paletteer_c("grDevices::SunsetDark", direction = -1, na.value = "white") +
  theme(text = element_text(size=9), axis.text.x = element_text(angle = -45, hjust = 0.1, vjust = 0.1),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=7),
        strip.background = element_rect(fill = "deeppink3", colour = "white"),
        strip.text = element_text(color = "white", size = 8),
        axis.text=element_text(size=7),
        axis.title=element_text(size = 8),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        plot.title = element_text(size = 10)) +
  facet_grid(~fB)

#pdf("FigureS3.pdf", width = 6.5, height = 6)
jpeg("FigureS3.jpg", width = 6.5, height = 6, units = "in", res = 300)
grid.arrange(g1, g2, nrow = 2)
dev.off()


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

SensResults_new <- SensAnalyze_short(fB_range, omega_b_range)
write.csv(SensResults_new, "SensResults_new.csv")

############## ANALYZE RESULTS ###############
SensResults_new <- read.csv("SensResults_new.csv")

sensres <- SensResults_new %>%
  filter(Profile != "S" & Profile != "prev" & type != "Inc" & type != "CumCases" & !is.na(Profile)) %>%
  group_by(time, ResistState, Strategy, fB, omega_b) %>%
  summarise(prev = sum(individuals)/10^6)

sensres <- sensres %>%
  ungroup() %>%
  group_by(Strategy, fB, omega_b, time) %>%
  mutate(prop = prev/sum(prev)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Strategy, fB, omega_b, time),
              names_from=ResistState,
              values_from = c(prev, prop)) 

colnames(sensres) <- c("Strategy", "fB", "omegaB", "time", "prev_AB", "prev_Aonly", "prev_Bonly", "prev_neither",
                       "prop_AB", "prop_Aonly", "prop_Bonly", "prop_neither")

sensres <- sensres %>%
  mutate(prev_all = prev_AB + prev_Aonly + prev_Bonly + prev_neither,
         prop_all = prop_AB + prop_Aonly + prop_Bonly + prop_neither,
         prop_res = prop_AB + prop_Aonly + prop_Bonly,
         prop_A = prop_AB + prop_Aonly,
         prop_B = prop_AB + prop_Bonly)

incidence <- SensResults_new %>%
  filter(type == "CumCases") %>%
  pivot_wider(id_cols = c(Strategy, fB, omega_b, time),
              names_from = type,
              values_from = individuals)
colnames(incidence)[colnames(incidence) == "omega_b"] <- "omegaB"

sensres <- left_join(sensres, incidence)

sensres <- sensres %>%
  mutate(Strategy_cat = rep(NA, nrow(.)))
sensres$Strategy_cat[sensres$Strategy == 1] <- "1. Random 50-50\nAllocation"
sensres$Strategy_cat[sensres$Strategy == 2] <- "2. Combination \nTherapy"
sensres$Strategy_cat[sensres$Strategy == 3] <- "3. Reserve \nStrategy"
sensres$Strategy_cat[sensres$Strategy == 4] <- "4. Gradual\nSwitch"

sensres <- sensres %>%
  group_by(fB, omegaB) %>%
  dplyr::mutate(ID = cur_group_id()) %>%
  ungroup()

sensres <- sensres %>%
  mutate(prop_B = ifelse(prop_B < 10^-6, 9e-7, prop_B))


####### NEW VIZ IDEA ########

year10 <- subset(sensres, time == 3650)

library(ggridges)
library(paletteer)

f1 <- ggplot(year10, aes(x = as.character(Strategy_cat), y = log10(prop_A))) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.4) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.title.y=element_blank()) +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = "Probability of de novo\nResistance, Drug B") +
  ylab("Proportion") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-4, -2, 0), labels = c("1e-04", 0.01, 1),
                     limits = c(-4.5, 0)) +
  coord_flip() +
  ggtitle("A. Year 10: Proportion of Infections\n     Resistant to Drug A")


f2 <- ggplot(year10, aes(x = as.character(Strategy), y = log10(prop_B))) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.3) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y=element_blank(),
        legend.position = "none") +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = "Probability of de novo\nResistance, Drug B") +
  ylab("Proportion") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-6, -4, -2, 0), labels = c("1e-06", "1e-04", 0.01, 1),
                     limits = c(-6.2, 0)) +
  coord_flip() +
  ggtitle("B. Year 10: Proportion of Infections\n     Resistant to Drug B")

f3 <- ggplot(year10, aes(x = as.character(Strategy), y = CumCases)) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.4) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        plot.title = element_text(face = "bold")) +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = "Probability of de novo\nResistance, Drug B") +
  ylab("Case Count") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(1780000, 1800000, 1820000, 1840000, 1860000), 
                     labels = c("1.78 M", "1.80 M", "1.82 M", "1.84 M", "1.86 M")) +
  coord_flip() +
  ggtitle("C. Year 10: Cumulative Number\n     of New Infections") 

###
year20 <- subset(sensres, time == 7300)

f4 <- ggplot(year20, aes(x = as.character(Strategy_cat), y = log10(prop_A))) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.3) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none",
        axis.title.y=element_blank()) +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = "Probability of de novo\nResistance, Drug B") +
  ylab("Proportion") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-4, -2, 0), labels = c("1e-04", 0.01, 1),
                     limits = c(-4.5, 0)) +
  coord_flip() +
  ggtitle("D. Year 20: Proportion of Infections\n     Resistant to Drug A")


f5 <- ggplot(year20,aes(x = as.character(Strategy), y = log10(prop_B))) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.3) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y=element_blank(),
        legend.position = "none") +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = "Probability of de novo\nResistance, Drug B") +
  ylab("Proportion") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-6, -4, -2, 0), labels = c("1e-06", "1e-04", 0.01, 1),
                     limits = c(-6.2, 0)) +
  coord_flip() +
  ggtitle("E. Year 20: Proportion of Infections\n     Resistant to Drug B")

f6 <- ggplot(year20, aes(x = as.character(Strategy), y = CumCases)) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.4) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        plot.title = element_text(face = "bold")) +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = "Probability of de novo\nResistance, Drug B") +
  ylab("Case Count") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(4000000, 5000000, 6000000), labels = c("4 M", "5 M", "6 M")) +
  coord_flip() +
  ggtitle("F. Year 20: Cumulative Number\n     of New Infections")

#pdf("SensResShortterm_new.pdf", width = 12.5, height = 8.5)
jpeg("SensResShortterm_new.jpg", width = 12.5, height = 8.5, units = "in", res = 300)
grid.arrange(f1, f2, f3, f4, f5, f6, nrow = 2, widths= c(1, 0.8, 1.2))
dev.off()
