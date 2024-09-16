library(tidyverse)
library(patchwork)
library(ggrepel)
library(geosphere)

#read in FOI file
#FOI for 365 days for an exemplar month

df.conc <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/foi_decay_exemplar.csv")


#Specifying some half lives for illustrating decay in protection
r_fast <- 0.010
r_mid  <- 0.005
r_slow <- 0.003
#Starting at 1
P_0 <- 1
#Decay of protection: Value gives the relative efficacy of the prevention measure for that date
#E.g., when protection = 0   then probability of infection = 1*prob.inf
#E.g., when protection = 0.5 then probability of infection = 0.5*prob.inf
#E.g., when protection = 1   then probability of infection = 0*prob.inf

#Scenarios:
#No supplemental intervention:  none
#Short half life, fast decay:   fast
#Mid half life:                 mid
#Long half life, slow decay:    slow
#Continuous, no decay:          cont # Here using estimates of 88% efficacy for SMC from:
# https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)32227-3/fulltext

df.decay <- tibble(day  = 1:(365*2),
                   #"norm" scenario shows growth in the absence of any protection
                   none = rep(1, 365*2),
                   fast = 1-P_0*(1 - r_fast)^day,
                   mid  = 1-P_0*(1 - r_mid)^day,
                   slow = 1-P_0*(1 - r_slow)^day,
                   cont = rep(0.12, 365*2))

#Note that mod = 1- protection (a modifier multiplied by the expected proportion infected)

df.conc2 <- left_join(df.conc, df.decay, by = join_by(day)) %>%
  pivot_longer(!c(day:prob.inf), names_to = "decay_scenario", values_to = "mod") %>%
  mutate(prob.inf.c = prob.inf*mod) %>%
  mutate(decay_scenario = factor(decay_scenario, levels = c("none", "fast", "mid", "slow", "cont"))) %>%
  #filter(site_code %in% c("S1", "S2", "N5")) %>%
  filter(site_code %in% c("S2", "N3")) %>%
  mutate(foi_cat = case_when(
    site_code == "S2" ~ "High",
    site_code == "N3" ~ "Low")) %>%
  mutate(foi_cat = factor(foi_cat, levels = c("High", "Low"))) %>%
  dplyr::select(-age) %>%
  mutate(interv_scenario = case_when(
    decay_scenario == "none" ~ "No supplemental intervention",
    decay_scenario == "fast" ~ "Shorter half-life intervention",
    decay_scenario == "mid"  ~ "Mid half-life intervention",
    decay_scenario == "slow" ~ "Longer half-life intervention",
    decay_scenario == "cont" ~ "Continuous intervention coverage")) %>%
  mutate(interv_scenario = factor(interv_scenario, levels = c("No supplemental intervention",
                                                              "Shorter half-life intervention",
                                                              "Mid half-life intervention",
                                                              "Longer half-life intervention",
                                                              "Continuous intervention coverage")))


#Threshhold crossing values
#Where each line hits a chosen threshhold (using 0.35 as an example)
v.threshhold <- 0.35

df.conc2.threshholds <- df.conc2 %>% 
  filter(decay_scenario != "cont") %>% 
  group_by(decay_scenario, foi_cat) %>% 
  filter(abs(prob.inf.c - v.threshhold) == min(abs(prob.inf.c - v.threshhold))) %>%
  mutate(threshold = v.threshhold)

#Plotting
p.curve <- df.conc2 %>% 
  filter(day < 500) %>%
  ggplot(aes(x = day, y = prob.inf.c, color = interv_scenario)) +
  geom_hline(yintercept = 0.35, color = "grey50", 
             #linetype = "dotted",
             alpha = 0.5,
             linewidth = 0.5) +
  geom_line(aes(linetype = foi_cat), linewidth = 0.8) +
  geom_point(data = df.conc2.threshholds,
             aes(x = day, y = threshold),
             color = "grey50") +
  #geom_vline(xintercept = df.conc2.threshholds$day) + 
  # geom_segment(data = df.conc2.threshholds,
  #              aes(x = day, xend = day, y = 0, yend = threshold),
  #              color = "grey50", 
  #              alpha = 0.5,
  #              linewidth = 0.5) +
  scale_color_viridis_d(option = "inferno", name = "Intervention Scenario", direction = -1, begin = 0.25, end = 0.95) +
  xlab("Time") + ylab("Cumulative probability of a new malaria infection") +
  scale_y_continuous(breaks = c(0, 1), position = "left") +
  scale_linetype_discrete(name = "Force of infection\n(FOI)") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.position = "right",
        legend.justification = "top")
p.curve





