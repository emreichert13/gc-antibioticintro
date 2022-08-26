
# Visualize results from "GCtransmission.Rmd"
##Area plots of infections over time, colored by Resistance Profile
# under baseline model assumptions

source("GCtransmission.Rmd")

#Random 50-50 Treatment
infections1 <- random_sim_long %>%
  filter(Profile != "S") %>%
  group_by(time, ResistState) %>%
  summarise(prev = sum(individuals)/10^6) %>%
  mutate(percent = prev/sum(prev), Strategy = "1. 50-50 Treatment",
         LossOne = 18.7, LossBoth = 19.2, PrevEq = "9.03%")

#Combination therapy
infections2 <- combo_sim_long %>%
  filter(Profile != "S") %>%
  group_by(time, ResistState) %>%
  summarise(prev = sum(individuals)/10^6) %>%
  mutate(percent = prev/sum(prev), Strategy = "2. Combination Therapy",
         LossOne = 19.9, LossBoth = 19.9, PrevEq = "9.03%")

#reserve strategy
infections3 <- reserve_sim_long %>%
  filter(Profile != "S") %>%
  group_by(time, ResistState) %>%
  summarise(prev = sum(individuals)/10^6) %>%
  mutate(percent = prev/sum(prev), Strategy = "3. Reserve Strategy",
         LossOne = 6.5, LossBoth = 13.9, PrevEq = "9.56%")

#gradual strategy
infections4 <- gradual_sim_long %>%
  filter(Profile != "S") %>%
  group_by(time, ResistState) %>%
  summarise(prev = sum(individuals)/10^6) %>%
  mutate(percent = prev/sum(prev), Strategy = "4. Gradual Switch",
         LossOne = 13.6, LossBoth = 18.2, PrevEq = "9.03%")

#Monotherapy, Drug A
infections5 <- monotherapy_sim_long_A %>%
  filter(Profile != "S") %>%
  group_by(time, ResistState) %>%
  summarise(prev = sum(individuals)/10^6) %>%
  mutate(percent = prev/sum(prev), Strategy = "Monotherapy, Drug A",
         LossOne = 6.5, LossBoth = NA, PrevEq = "10.42%")

#Reserve, followed by combination therapy
infections6 <- reserve_combo_sim_long %>%
  filter(Profile != "S") %>%
  group_by(time, ResistState) %>%
  summarise(prev = sum(individuals)/10^6) %>%
  mutate(percent = prev/sum(prev), Strategy = "Reserve Strategy,\n followed by Combination",
         LossOne = 6.5, LossBoth = 13.9, PrevEq = "9.04%")

#join all results
infections <- rbind(infections1, infections2, infections3, infections4, infections5, infections6)

#plot prevalence of gonorrhea over time, by Resistance Profile
ggplot(infections,aes(x=time/365, y=prev*100, fill = ResistState)) +
  geom_area() + scale_fill_manual(values = c("deeppink3","salmon", "lightpink", "midnightblue")) +
  theme_classic() + facet_wrap(~Strategy) +
  xlab("Years") + ylab("Prevalence") + theme(text = element_text(size=19)) +
  geom_vline(aes(xintercept = LossOne), col = "white", linetype = "dashed") + 
  geom_vline(aes(xintercept = LossBoth), col = "white") + xlim(c(0,40))

# % of gonorrhea infections made up of each Resistance Profile, over time
pdf(file="AreaPlot.pdf", width=10, height=7)
#jpeg(file="AreaPlot.jpg", width=10, height=7, units = "in", res = 300)
ggplot(infections,aes(x=time/365, y=percent*100, fill = factor(ResistState, levels = c("Neither", "A only", "B only", "A and B")))) +
  geom_area() + scale_fill_manual(values = c("midnightblue","salmon", "lightpink", "deeppink3")) +
  theme_classic() + facet_wrap(~Strategy, scales = "free_x") +
  xlab("Years") + ylab("% of Gonorrhea Infections") + theme(text = element_text(size=15))  +
  geom_segment(aes(x = LossOne, xend = LossOne, y = 0, yend = 100), col = "salmon")+
  geom_segment(aes(x = LossBoth, xend = LossBoth, y = 0, yend = 100), col = "white")+
  xlim(c(0,40)) +
  geom_hline(aes(yintercept = 5), col = "white", linetype= "dashed") +
  labs(fill = "Resistance Profile") +
  geom_label(data = infections %>% filter(time == 1), aes(label = PrevEq), x = 33.5, y = 90,
             fill = "white", col = "black", label.size = 0.08)
dev.off()
