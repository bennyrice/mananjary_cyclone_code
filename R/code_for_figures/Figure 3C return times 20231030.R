library(tidyverse)
library(patchwork)
library(ggrepel)


#Reading in data (corresponding to function call in determining FOI script)
df.foi.i <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_for_return_time_plotting.csv")


#Only looking at school aged children
df.foi <- df.foi.i %>% filter(age == 13)

#Dropping AND bumping FOI by fold changes
df.foi.mod <- df.foi %>%
  mutate(foi.up.100   = prob.inf*2.00,
         foi.up.90    = prob.inf*1.90,
         foi.up.80    = prob.inf*1.80,
         foi.up.70    = prob.inf*1.70,
         foi.up.60    = prob.inf*1.60,
         foi.up.50    = prob.inf*1.50,
         foi.up.40    = prob.inf*1.40,
         foi.up.30    = prob.inf*1.30,
         foi.up.20    = prob.inf*1.20,
         foi.up.10    = prob.inf*1.10,
         foi.0        = prob.inf,      #baseline
         foi.down.10  = prob.inf*0.90,
         foi.down.20  = prob.inf*0.80,
         foi.down.30  = prob.inf*0.70,
         foi.down.40  = prob.inf*0.60,
         foi.down.50  = prob.inf*0.50,
         foi.down.60  = prob.inf*0.40,
         foi.down.70  = prob.inf*0.30,
         foi.down.80  = prob.inf*0.20,
         foi.down.90  = prob.inf*0.10,
         foi.down.99  = prob.inf*0.01)

df.foi.mod <- df.foi.mod %>%
  pivot_longer(!c(day:code.new), names_to = "foi.scenario", values_to = "foi.c") %>%
  dplyr::select(site_code, code.new, age, foi.scenario, foi.c)


f.return.time.i <- function(foi, threshold){
  
  v.prob.inf <- rep(NA, 10000)
  for(t in 1:10000){
    v.prob.inf[t] <- 1-exp(-foi*t)
  }
  return.time <- length(which(v.prob.inf < threshold))
  #If less than 1 day, will give zero which complicates plotting on log scale
  #If RT is between day 0 and day 1 then representing as day 1
  return.time <- ifelse(return.time < 1, 1, return.time)
  
  return(return.time)
}

df.foi.mod.RT <- df.foi.mod %>% rowwise() %>% mutate(RT  = f.return.time.i(foi = foi.c, threshold = 0.10))

df.foi.mod.RT.p <- df.foi.mod.RT %>% mutate(foi.scenario = factor(foi.scenario, levels = unique(df.foi.mod.RT$foi.scenario))) %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10"))) %>%
  mutate(age_cat = case_when(age >=  6 & age < 14 ~ "School aged children (13y)")) %>%
  mutate(foi.scenario = case_when(
    foi.scenario == "foi.up.100"  ~ "+100%",
    foi.scenario == "foi.up.90"   ~ "+90%",
    foi.scenario == "foi.up.80"   ~ "+80%",
    foi.scenario == "foi.up.70"   ~ "+70%",
    foi.scenario == "foi.up.60"   ~ "+60%",
    foi.scenario == "foi.up.50"   ~ "+50%",
    foi.scenario == "foi.up.40"   ~ "+40%",
    foi.scenario == "foi.up.30"   ~ "+30%",
    foi.scenario == "foi.up.20"   ~ "+20%",
    foi.scenario == "foi.up.10"   ~ "+10%",
    foi.scenario == "foi.0"       ~ "0%",
    foi.scenario == "foi.down.10" ~ "-10%",
    foi.scenario == "foi.down.20" ~ "-20%",
    foi.scenario == "foi.down.30" ~ "-30%",
    foi.scenario == "foi.down.40" ~ "-40%",
    foi.scenario == "foi.down.50" ~ "-50%",
    foi.scenario == "foi.down.60" ~ "-60%",
    foi.scenario == "foi.down.70" ~ "-70%",
    foi.scenario == "foi.down.80" ~ "-80%",
    foi.scenario == "foi.down.90" ~ "-90%",
    foi.scenario == "foi.down.99" ~ "-99%")) %>%
  mutate(foi.scenario = factor(foi.scenario, levels = c("+100%","+90%","+80%","+70%","+60%","+50%","+40%","+30%","+20%","+10%","0%","-10%","-20%","-30%","-40%","-50%","-60%","-70%","-80%","-90%","-99%")))

##########################################################################################
## Figure 3C: Plotting FOI scenarios (x) vs Days til infection reaches 10% (y)
##########################################################################################

df.foi.mod.RT.p %>% 
  ggplot(aes(x = foi.scenario, y = RT, group = code.new, color = code.new)) +
  #Shading portion of the plot to higlight scenarios < 90 days
  geom_rect(xmin=0.5, xmax=21.5, ymin=log10(3), ymax=log10(90), fill = "white", alpha = 0, linewidth = 1, linetype = 2, colour = "#EAC0BD") +
  geom_point(size = 1) +
  geom_path() +
  geom_hline(yintercept = 90, color = "grey50", alpha = 0.6, linetype = "11") +
  geom_vline(xintercept = "0%", color = "grey50", alpha = 0.6, linetype = "11") +
  xlab("% change relative to observed FOI\nForce of infection (FOI) scenario") + 
  ylab("Number of days until infection prevalence reaches threshold (10%)") +
  scale_color_viridis_d(option = "turbo", name = "Site Code") +
  facet_wrap(vars(age_cat), nrow = 1) +
  scale_y_continuous(trans = "log10",
                     breaks = c(1, 5, 10, 100, 1000, 10000)) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, family = "Andale Mono", hjust = 1, vjust = 0.5, size = 11),
        axis.title = element_text(family = "Arial", size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(family = "Arial", size = 14),
        panel.grid = element_blank())


