
library(ggpubr)
library(sigmoid)
library(gridExtra)

# This script is used to generate Supplementary Fig.1 for Resistance-minimizing strategies for introducing 
# a novel antibiotic for gonorrhea treatment: a mathematical modeling study by Reichert et al.

#Reserve Strategy
 # explore transition between drug A and drug B once drug A hits 5% resistance threshold

# Function will transition between 0 and 1 when t and Q are approximately equal
smooth.transition1_new <- function(t, Q, tune = 0.02){
  sigmoid((t/Q -1)/tune)
}
smooth.transition2_new <- function(t, Q, tune = 0.02){
  1-sigmoid((t/Q - 1)/tune)
}

t <- c(seq(0, 365*40, 1))
Q <- 365*10 + 182.5 # if drug A hits 5% at 10 years, for illustration
prA_1 <- smooth.transition2_new(t, Q)
prB_1 <- smooth.transition1_new(t, Q)

p1 <- ggplot() + 
  geom_point(aes(x = t/365, y = prA_1), col = "navy") + 
  geom_point(aes(x = t/365, y = prB_1), col = "deeppink") +
  theme_classic() + 
  ylab("Probability of Treatment with Drug") +
  xlab("Time (Years)") +
  annotate("text", x = 35, y = 0.95, label = "B") +
  annotate("text", x = 35, y = 0.05, label = "A") + 
  geom_vline(xintercept = 10, linetype = "dashed") +
  scale_x_continuous(breaks = c(10, 40), labels = c("5% Threshold Met", 40)) +
  ggtitle("A. Reserve strategy")

#Gradual Switch
# Sigmoid transition with midpoint of 4 years
# Function will transition between 0 and 0.5 with midpoint of Q

t <- c(seq(0, 365*40, 1))
smooth.transition3 <- function(prevA, Q, tune = .22){
  sigmoid((t/Q-1)/tune)/2
}
smooth.transition4 <- function(prevA, Q, tune = .22){
  1-sigmoid((t/Q-1)/tune)/2
}
Q <- 365*4
prA_2 <- smooth.transition4(t, Q)
prB_2 <- smooth.transition3(t, Q)


p2 <- ggplot() + 
  geom_point(aes(x = t/365, y = prA_2), col = "navy") + 
  geom_point(aes(x = t/365, y = prB_2), col = "deeppink") + theme_classic() + 
  xlab("Time (Years)") +
  ylab("Probability of Treatment with Drug") +
  annotate("text", x = 5, y = 0.1, label = "B") +
  annotate("text", x = 5, y = 0.9, label = "A") + 
  ggtitle("B. Gradual Switch strategy")

prAB_T <- data.frame(cbind(t = t, prA = round(prA_2, 3), prB = round(prB_2,3)))

#jpeg(file="Supp_transitions_update.jpg", width=8, height=4, units = "in", res = 300)
#pdf("Supp_transitions_update.pdf", width = 8, height = 4)
figure <- grid.arrange(p1, p2, nrow = 1)
figure
#dev.off()

# #Gradual Switch, can play around with faster or slower versions
# t <- c(seq(0, 365*40, 1))
# smooth.transition3 <- function(prevA, Q, tune = .3){
#   sigmoid((t/Q-1)/tune)/2
# }
# smooth.transition4 <- function(prevA, Q, tune = .3){
#   1-sigmoid((t/Q-1)/tune)/2
# }
# Q <- 365*2
# prA_3 <- smooth.transition4(t, Q)
# prB_3 <- smooth.transition3(t, Q)
# 
# p3 <- ggplot() + 
#   geom_point(aes(x = t/365, y = prA_3), col = "navy") + 
#   geom_point(aes(x = t/365, y = prB_3), col = "deeppink") + theme_classic() + 
#   xlab("Time (Years)") +
#   ylab("Probability of Treatment with Drug") +
#   annotate("text", x = 4, y = 0.1, label = "B") +
#   annotate("text", x = 4, y = 0.9, label = "A") + 
#   ggtitle("C. Gradual introduction (fast)")

#4 year midpoint - equilibrium reached ~ 9.4 years
#3 year midpoint - equilibrium reached ~ 7.0 years
#2 year midpoint - equilibrium reached ~ 4.7 years
#1 year midpoint - equilibrium reached ~ 2.4 years




