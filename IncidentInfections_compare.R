source("GCtransmission.Rmd")

# This script calculates the number of incident infections each day, under each strategy
# as well as cumulative case counts 

# Then, makes visual comparisons in cumulative cases averted and IRRs by strategy
# under baseline model parameters

# I. Make our dataset

# Strategy 1 - Random 50-50 Allocation
parms$prA <- 0.5          #Pr of treatment with A
parms$prB <- 0.5          #Pr of treatment with B

random_inc <- as.data.frame(ode(inits, dt, random.SI, parms=parms)) %>%
  mutate(S = S1 + S2 + S3,
         Strategy = "1. 50-50 Treatment",
         CumCases = cumsum(Inc)) %>%
  select(time, S, Inc, CumCases, Strategy)

# Monotherapy, for comparison
parms$prA <- 1          #Pr of treatment with A
parms$prB <- 0          #Pr of treatment with B

monotherapy_inc <- as.data.frame(ode(inits, dt, random.SI, parms=parms)) %>%
  mutate(S = S1 + S2 + S3,
         Strategy = "Monotherapy, Drug A",
         CumCases = cumsum(Inc)) %>%
  select(time, S, Inc, CumCases, Strategy)

# Strategy 2 - Combination Therapy
combo_inc <- as.data.frame(ode(inits, dt, combo.SI, parms=parms)) %>%
  mutate(S = S1 + S2 + S3,
         Strategy = "2. Combination Therapy",
         CumCases = cumsum(Inc)) %>%
  select(time, S, Inc, CumCases, Strategy)  

# Strategy 3 - Reserve Strategy
reserve_inc <- as.data.frame(ode(inits, dt, reserve.SI, parms=parms)) %>%
  mutate(S = S1 + S2 + S3,
         Strategy = "3. Reserve Strategy",
         CumCases = cumsum(Inc)) %>%
  select(time, S, Inc, CumCases, Strategy) 

# Strategy 4 - Gradual Switch
gradual_inc <- as.data.frame(ode(inits, dt, gradual.SI, parms=parms)) %>%
  mutate(S = S1 + S2 + S3,
         Strategy = "4. Gradual Switch",
         CumCases = cumsum(Inc)) %>%
  select(time, S, Inc, CumCases, Strategy) 

# Reserve followed by Combination Therapy, for comparison
reserve_combo_inc <- as.data.frame(ode(inits, dt, reserve_combo.SI, parms=parms, method = "ode45")) %>%
  mutate(S = S1 + S2 + S3,
         Strategy = "Reserve Strategy,\n followed by Combination",
         CumCases = cumsum(Inc)) %>%
  select(time, S, Inc, CumCases, Strategy)

# join all results
incidence <- rbind(random_inc, monotherapy_inc, combo_inc, reserve_inc, gradual_inc, reserve_combo_inc)

# II. Analyze results

incidence %>% group_by(Strategy) %>% 
  summarise(max(CumCases))

incidence %>% group_by(Strategy) %>%
  summarise(sum(Inc))

#calculate avg. no of incident infections per year, under these strategies
#throughout lifespan of the drugs, A and/or B
summ_cases <- incidence %>%
  mutate(Year = time/365,
         LossBoth = ifelse(Strategy == "1. 50-50 Treatment", 19.2,
                           ifelse(Strategy == "2. Combination Therapy", 19.9,
                                  ifelse(Strategy == "3. Reserve Strategy", 13.9,
                                         ifelse(Strategy == "4. Gradual Switch", 18.2, 
                                                ifelse(Strategy == "Monotherapy, Drug A", 6.5, NA)))))) %>%
  filter(Year <= LossBoth) %>%
  group_by(Strategy, LossBoth) %>%
  summarise(NewCases = sum(Inc), NewCases2 = max(CumCases)) %>%
  mutate(AvgYr = NewCases/LossBoth)

#Rolling incidence per year, to look at max and min
cases_yr <- incidence %>%
  mutate(Year = time/365,
         LossBoth = ifelse(Strategy == "1. 50-50 Treatment", 19.2,
                           ifelse(Strategy == "2. Combination Therapy", 19.9,
                                  ifelse(Strategy == "3. Reserve Strategy", 13.9,
                                         ifelse(Strategy == "4. Gradual Switch", 18.2, 
                                                ifelse(Strategy == "Monotherapy, Drug A", 6.5, NA)))))) %>%
  filter(Year <= LossBoth) %>%
  group_by(Strategy) %>%
  arrange(time) %>%
  mutate(Infections_RollYr = CumCases - lag(CumCases, 365)) %>%
  summarise(min(Infections_RollYr, na.rm = T), 
            max(Infections_RollYr, na.rm = T))

#Calculate Incidence Rate, per 100 PY
incidence <- incidence %>%
  mutate(IR = Inc/S * 365 * 100)

compare <- incidence %>%
  filter(Strategy == "3. Reserve Strategy") %>%
  select(time, IR_compare = IR)

#IRRs and IRDs, relative to reserve strategy
incidence <- left_join(incidence, compare) %>%
  mutate(IRR = IR/IR_compare,
         IRD = IR - IR_compare)

compare <- incidence %>%
  filter(Strategy == "3. Reserve Strategy") %>%
  select(time, CumCases_compare = CumCases)

#Cases Averted, relative to reserve strategy
incidence <- left_join(incidence, compare) %>%
  mutate(CasesAverted = (CumCases_compare - CumCases))

incidence <- incidence %>%
  mutate(LossBoth = ifelse(Strategy == "1. 50-50 Treatment", 19.2,
                           ifelse(Strategy == "2. Combination Therapy", 19.9,
                                  ifelse(Strategy == "3. Reserve Strategy", 13.9,
                                         ifelse(Strategy == "4. Gradual Switch", 18.2, 
                                                ifelse(Strategy == "Monotherapy, Drug A", 6.5, NA))))))

incidence <- incidence %>%
  mutate(time = time/365) %>%
  filter(Strategy != "Monotherapy, Drug A") 

# Plot IRRs over time under baseline model
f1 <- ggplot(data = incidence) +
  geom_line(aes(x = time, y = IRR, col = Strategy), size = 2) +
  ylab("Incidence Rate Ratio") +
  geom_vline(aes(xintercept = LossBoth, col = Strategy), linetype = "dashed", size = 0.8) +
  scale_discrete_manual("color", values = 
                          paletteer_d("ggthemes::excel_Ion_Boardroom")) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Years") +
  ggtitle("A.")

# Plot cumulative cases averted (relative to Reserve strategy) over time, under baseline model
f2 <- ggplot(data = incidence) +
  geom_line(aes(x = time, y = CasesAverted, col = Strategy), size = 2) +
  ylab("Cumulative Cases Averted, \n relative to Reserve Strategy") +
  geom_vline(aes(xintercept = LossBoth, col = Strategy), linetype = "dashed", size = 0.8)  +
  scale_discrete_manual("color", values = 
                          paletteer_d("ggthemes::excel_Ion_Boardroom")) +
  scale_y_continuous(breaks = c(500000, 1000000, 1500000), labels = c("0.5 M", "1.0 M", "1.5 M")) +
  theme_classic() +
  xlab("Years") +
  ggtitle("B.")

#export plot
pdf("IncidenceFig.pdf", width = 11, height = 4)
#jpeg("IncidenceFig.jpg", width = 11, height = 4, units = "in", res = 300)
grid.arrange(f1, f2, nrow = 1, widths = c(0.4, 0.6))
dev.off()
  