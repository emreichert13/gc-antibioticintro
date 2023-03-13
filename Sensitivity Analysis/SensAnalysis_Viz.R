
#E Reichert 2022
#Analyze/visualize results of sensitivity analysis produced by "SensAnalysis_Run.R"

## Used to generate Figs. 2-3 and Supp. Fig 4 for Resistance-minimizing strategies for introducing 
# a novel antibiotic for gonorrhea treatment: a mathematical modeling study by Reichert et al.

library(ggridges)
library(paletteer)
library(ggplot2)
library(readr)
library(dplyr)

################ ANALYZE LONG TERM RESULTS ####################

## I - Baseline Model Scenario
sens_baseline <- read.csv("SensResults_Baseline.csv") %>% select(-X)

SensResults <- sens_baseline

SensResults$LossA[is.infinite(SensResults$LossA) | SensResults$LossA > 100] <- 100
SensResults$LossB[is.infinite(SensResults$LossB)| SensResults$LossB > 100] <- 100
SensResults$LossAB[is.infinite(SensResults$LossAB) | SensResults$LossAB > 100] <- 100
SensResults <- SensResults %>%
  filter(Strategy != 0 & Strategy != 4.2)

SensResults$LossBoth <- apply(SensResults[, 4:5], 1, max)
SensResults$LossOne <- apply(SensResults[, 4:5], 1, min)

SensResults$Strategy[SensResults$Strategy == 1] <- "1. 50-50"
SensResults$Strategy[SensResults$Strategy == 2] <- "2. Combination"
SensResults$Strategy[SensResults$Strategy == 3] <- "3. Reserve"
SensResults$Strategy[SensResults$Strategy == 4.1] <- "4. Gradual"

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
  mutate(LossBoth = ifelse(fB == 0.80 & omega_B <= 10^-6, NA, LossBoth),
         Diff_both = ifelse(fB == 0.80 & omega_B <= 10^-6, NA, Diff_both))

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

my_y_title <- expression(paste("Probability of ", italic("de novo"), " Resistance, Drug B"))
#pdf("Figure2.pdf", width = 7.5, height = 3.5)
#jpeg("Figure2.jpg", width = 7.5, height = 3.5, units = "in", res = 300)
ggplot(SensResults, aes(Strategy, reorder(omega_B, desc(omega_B)), fill= Diff_both)) + 
  geom_tile() + theme_classic() +
  geom_text(aes(label = Label), col = "white", size = 2.5) +
  xlab("Strategy") + 
  ylab(my_y_title) + 
  labs(fill = "Additional Time to Loss \nof Both Drugs, relative to \nReserve Strategy (years)") + 
  scale_fill_paletteer_c("grDevices::SunsetDark", direction = -1, na.value = "#660066") +
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
#dev.off()

## II - Reduced fitness associated with drug A resistance (fA = 0.90)

sens_fA <- read.csv("SensResults_fA.csv") %>% select(-X)

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

sens_fA$Label[sens_fA$Label == ">NA"] <- NA

g1 <- ggplot(sens_fA, aes(Strategy, reorder(omega_B, desc(omega_B)), fill= Diff_both)) + 
  geom_tile() + theme_classic() +
  geom_text(aes(label = Label), col = "white", size = 2.5) +
  xlab("Strategy") + 
  ylab(my_y_title) + 
  labs(fill = "Additional Time to Loss \nof Both Drugs, relative to \nReserve Strategy (years)",
       title = "A. Relative Fitness of Drug A Resistant Strains (fA) = 0.90",
       subtitle = "fB = Relative Fitness of Drug B Resistant Strains") + 
  scale_fill_paletteer_c("grDevices::SunsetDark", direction = -1, na.value = "#660066") +
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

## III - Higher Pr of de novo resistance associated with Drug A treatment (omegaA = 10^-4)

sens_omegaA <- read.csv("SensResults_omegaA.csv") %>% select(-X)

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

sens_omegaA$Label[sens_omegaA$Label == ">NA"] <- NA

g2 <- ggplot(sens_omegaA, aes(Strategy, reorder(omega_B, desc(omega_B)), fill= Diff_both)) + 
  geom_tile() + theme_classic() +
  geom_text(aes(label = Label), col = "white", size = 2.5) +
  xlab("Strategy") + 
  ylab(my_y_title) + 
  labs(fill = "Additional Time to Loss \nof Both Drugs, relative to \nReserve Strategy (years)",
       title = expression(paste("B. Probability of ", italic("de novo"), " Resistance, Drug A = 1e-04")),
       subtitle = "fB = Relative Fitness of Drug B Resistant Strains") + 
  scale_fill_paletteer_c("grDevices::SunsetDark", direction = -1, na.value = "#660066") +
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

library(gridExtra)
#pdf("FigureS3.pdf", width = 6.5, height = 6)
#jpeg("FigureS3.jpg", width = 6.5, height = 6, units = "in", res = 300)
grid.arrange(g1, g2, nrow = 2)
#dev.off()


############## ANALYZE SHORT TERM RESULTS ###############
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

# Visualize distribution of outcomes at 10 years
year10 <- subset(sensres, time == 3650)

summary(year10$prop_A)
summary(year10$prop_B)
summary(year10$CumCases)

year10 %>%
  group_by(Strategy_cat) %>%
  summarise(median(CumCases))

my_axis_title <- expression(paste("Probability of Resistance\nupon Treatment, Drug B"))

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
                        name = my_axis_title) +
  ylab("Proportion (log10 scale)") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-4, -2, 0), labels = c("1e-04", 0.01, 1), limits = c(-5, 0.1)) +
  coord_flip() +
  ggtitle("A. Year 10: Infections\n     Resistant to Drug A")


f2 <- ggplot(year10, aes(x = as.character(Strategy), y = log10(prop_B))) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.4) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y=element_blank(),
        legend.position = "none") +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = my_axis_title) +
  ylab("Proportion (log10 scale)") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-6, -4, -2, 0), labels = c("1e-06", "1e-04", 0.01, 1), limits = c(-6.1, 0.1)) +
  coord_flip() +
  ggtitle("B. Year 10: Infections\n     Resistant to Drug B")

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
                        name = my_axis_title) +
  ylab("Case Count") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(1780000, 1820000, 1860000), 
                     labels = c("1.78 M", "1.82 M", "1.86 M")) +
  coord_flip() +
  ggtitle("C. Year 10: Cumulative \n     Infections") 

# Visualize distribution of outcomes at 20 years
year20 <- subset(sensres, time == 7300)

summary(year20$prop_A)
summary(year20$prop_B)
summary(year20$CumCases)

year20 %>%
  group_by(Strategy_cat) %>%
  summarise(median(CumCases), mean(CumCases))

f4 <- ggplot(year20, aes(x = as.character(Strategy_cat), y = log10(prop_A))) +
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
                        name = my_axis_title) +
  ylab("Proportion (log10 scale)") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-4, -2, 0), labels = c("1e-04", 0.01, 1), limits = c(-5, 0.1)) +
  coord_flip() +
  ggtitle("D. Year 20: Infections\n     Resistant to Drug A")


f5 <- ggplot(year20,aes(x = as.character(Strategy), y = log10(prop_B))) +
  geom_boxplot(fill = "wheat1", col = "sandybrown", outlier.color = "white", width = 0.4) +
  geom_jitter(width = 0.1, height = 0.05,
              shape = 16, size = 2.5,
              aes(color = factor(omegaB))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y=element_blank(),
        legend.position = "none") +
  scale_discrete_manual("color", values = 
                          paletteer_d("vapoRwave::newRetro"),
                        name = my_axis_title) +
  ylab("Proportion (log10 scale)") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(-6, -4, -2, 0), labels = c("1e-06", "1e-04", 0.01, 1), limits = c(-6.1, 0.1)) +
  coord_flip() +
  ggtitle("E. Year 20: Infections\n     Resistant to Drug B")

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
                        name = my_axis_title) +
  ylab("Case Count") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_continuous(breaks = c(4000000, 5000000, 6000000), labels = c("4 M", "5 M", "6 M")) +
  coord_flip() +
  ggtitle("F. Year 20: Cumulative \n     Infections")

#pdf("SensResShortterm_new.pdf", width = 11, height = 7.5)
#jpeg("SensResShortterm_new.jpg", width = 11, height = 7.5, units = "in", res = 300)
grid.arrange(f1, f2, f3, f4, f5, f6, nrow = 2, widths= c(1, 0.8, 1.2))
#dev.off()
