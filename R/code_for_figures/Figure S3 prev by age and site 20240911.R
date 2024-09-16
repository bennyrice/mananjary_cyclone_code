library(tidyverse)
library(patchwork)
library(ggrepel)


########################################################################################################
## Reading in data (pre-processed RDT file in long format)

df5 <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_malaria_rdt/rdt_data_20231007.csv")

########################################################################################################
#Confidence intervals: using binom.test() and Clopper and Pearson (1934) method

#Writing a function to run binom.test() using number of positives and sample size
#Output of binom.test() is a htest object
#From which we extract the upper and lower confidence interval bounds
#For lower limit first
f.CI_calc_lower <- function(v.n_pos, v.n){
  v.CI_lower <- rep(NA, length(v.n_pos))
  for (i in 1:length(v.n_pos)) {
    x.i <- v.n_pos[i]
    n.i <- v.n[i]
    htest.i <- binom.test(x.i, n.i)
    v.CI_lower[i] <- htest.i$conf.int[1]
  }
  return(v.CI_lower)
}
#And now for upper limit
f.CI_calc_upper <- function(v.n_pos, v.n){
  v.CI_upper <- rep(NA, length(v.n_pos))
  for (i in 1:length(v.n_pos)) {
    x.i <- v.n_pos[i]
    n.i <- v.n[i]
    htest.i <- binom.test(x.i, n.i)
    v.CI_upper[i] <- htest.i$conf.int[2]
  }
  return(v.CI_upper)
}

########################################################################################################
#Re-labeling time points to reflect baseline as T0
df5 <- df5 %>% mutate(time_point = case_when(
  time_point == "T01" ~ "T0",
  time_point == "T02" ~ "T01",
  time_point == "T03" ~ "T02",
  time_point == "T04" ~ "T03",
  time_point == "T05" ~ "T04",
  time_point == "T06" ~ "T05",
  time_point == "T07" ~ "T06",
  time_point == "T08" ~ "T07",
  time_point == "T09" ~ "T08",
  time_point == "T10" ~ "T09",
  time_point == "T11" ~ "T10"
))

########################################################################################################
#Adding udpated site code labels
df5 <- df5 %>% mutate(site_code = case_when(
  site_code == "S1" ~ "MNJ.06",
  site_code == "S2" ~ "MNJ.07",
  site_code == "S3" ~ "MNJ.08",
  site_code == "S5" ~ "MNJ.09",
  site_code == "S6" ~ "MNJ.10",
  site_code == "N1" ~ "MNJ.05",
  site_code == "N2" ~ "MNJ.04",
  site_code == "N3" ~ "MNJ.03",
  site_code == "N4" ~ "MNJ.02",
  site_code == "N5" ~ "MNJ.01")) %>%
  mutate(site_code = factor(site_code, levels = c("MNJ.01",
                                                  "MNJ.02",
                                                  "MNJ.03",
                                                  "MNJ.04",
                                                  "MNJ.05",
                                                  "MNJ.06",
                                                  "MNJ.07",
                                                  "MNJ.08",
                                                  "MNJ.09",
                                                  "MNJ.10")))
  



########################################################################################################

#Pulling out numbers for summary stats
df.sum.age <- df5 %>% dplyr::select(-enrollment_status) %>% drop_na() %>%
  filter(age.yrs.at.sample >= 6) %>% filter(age.yrs.at.sample < 14) %>%
  group_by(time_point, site_code, age.cat.at.sample) %>%
  summarize(n = n(), n_pos = sum(rdt.result), prev = n_pos/n*100) %>% ungroup()

########################################################################################################



########################################################################################################
## Figure S3A ## Prevalence at baseline
########################################################################################################


#First, all ages combined
#Dropping rows without RDT data
#baseline (T0)
df.b1 <- df5 %>% dplyr::select(-enrollment_status) %>% drop_na() %>%
  group_by(time_point, site_code) %>%
  summarize(n = n(), n_pos = sum(rdt.result), prev = n_pos/n*100) %>% ungroup() %>%
  filter(time_point == "T0")
#Adding confidence intervals
#Multiplying by 100 to make a percentage
df.b1 <- df.b1 %>% rowwise() %>% 
  mutate(Prev.CI.lower = f.CI_calc_lower(n_pos, n) * 100) %>%
  mutate(Prev.CI.upper = f.CI_calc_upper(n_pos, n) * 100)

#Plotting
p.b1 <- df.b1 %>% 
  ggplot(aes(x=time_point, y=prev, color=site_code)) + 
  geom_point(stat = "identity", size = 3) +
  #Adding error bars
  geom_errorbar(aes(ymin=Prev.CI.lower, ymax=Prev.CI.upper, color=site_code),
                linewidth = 0.3,     # Thinner lines
                width = 0.06,        # Thinner horizontals on the Ts
                position="identity") +
  scale_color_manual(values = site_plot_colors1) +
  labs(title = "Baseline malaria prevalence", subtitle = "All ages") +
  ylab("Percent of sampled individuals positive by RDT (%)") +
  xlab("Locality") +
  ylim(0, 100) +
  facet_wrap(vars(site_code), nrow = 1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(), #Dropping gridlines
        panel.grid.minor = element_blank(), #Dropping gridlines
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p.b1


#School aged children at baseline
df.b1.sac <- df5 %>% dplyr::select(-enrollment_status) %>% drop_na() %>%
  filter(time_point == "T0") %>%
  filter(age.cat.at.sample == "C_6_13ys") %>%
  group_by(time_point, site_code) %>%
  summarize(n = n(), n_pos = sum(rdt.result), prev = n_pos/n*100) %>% ungroup() %>% 
  rowwise() %>% 
  mutate(Prev.CI.lower = f.CI_calc_lower(n_pos, n) * 100) %>%
  mutate(Prev.CI.upper = f.CI_calc_upper(n_pos, n) * 100)
#Plotting
p.b1.sac <- df.b1.sac %>% 
  ggplot(aes(x=time_point, y=prev)) + 
  geom_point(stat = "identity", size = 3, color="#EE564B") +
  #Adding error bars
  geom_errorbar(aes(ymin=Prev.CI.lower, ymax=Prev.CI.upper),
                color="#EE564B",
                linewidth = 0.3,     # Thinner lines
                width = 0.06,        # Thinner horizontals on the Ts
                position="identity") +
  labs(title = "Baseline malaria prevalence", subtitle = "School aged children (6-13 years)") +
  ylab("Percent of sampled individuals positive by RDT (%)") +
  xlab("Locality") +
  ylim(0, 100) +
  facet_wrap(vars(site_code), nrow = 1) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(), #Dropping gridlines
        panel.grid.minor = element_blank(), #Dropping gridlines
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p.b1.sac

p.b1 + p.b1.sac

########################################################################################################
## Figure S3B ## Prevalence by site and time point (all ages)
########################################################################################################

#Plot 
site_plot_colors1 <- c(rep("firebrick3", 10))
#Dropping rows without RDT data
df.ts1 <- df5 %>% dplyr::select(-enrollment_status) %>% drop_na() %>%
  group_by(time_point, site_code) %>%
  summarize(n = n(), n_pos = sum(rdt.result), prev = n_pos/n*100) %>% ungroup()
#Adding confidence intervals
#Multiplying by 100 to make a percentage
df.ts1 <- df.ts1 %>% rowwise() %>% 
  mutate(Prev.CI.lower = f.CI_calc_lower(n_pos, n) * 100) %>%
  mutate(Prev.CI.upper = f.CI_calc_upper(n_pos, n) * 100)
#Plotting
#Dropping baseline
p.ts1b <- df.ts1 %>% 
  filter(time_point != "T0") %>%
  ggplot(aes(x=time_point, y=prev, color=site_code)) + 
  #Adding cyclone time points
  geom_vline(xintercept = c(2.5, 8.5), alpha = 0.5, color = "grey20", linetype = "dashed") +
  #geom_vline(xintercept = "T09", alpha = 0.5, color = "grey20", linetype = "dashed") +
  geom_point(stat = "identity") +
  #Adding error bars
  geom_errorbar(aes(ymin=Prev.CI.lower, ymax=Prev.CI.upper, color=site_code),
                linewidth = 0.2,     # Thinner lines
                width = 0.05,        # Thinner horizontals on the Ts
                position="identity") +
  scale_color_manual(values = site_plot_colors1) +
  xlab("Time Point") + 
  ylab("Percent of sampled individuals positive by RDT (%)") +
  labs(title = "Malaria prevalence by time point", subtitle = "All ages") +
  facet_wrap(vars(site_code), nrow = 2) +
  theme_bw() + 
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(), #Dropping gridlines
        panel.grid.minor = element_blank(), #Dropping gridlines
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, family = "Andale Mono"))
p.ts1b


########################################################################################################
## Figure S3c ## Age distribution before and after cyclones
########################################################################################################

#Batsirai
#Original numbering: T02-T03 vs T04-T05
#New numbering:      T01-T02 vs T03-T04
#Freddy
#Original numbering: T08-T09 vs T10-T11
#New numbering:      T07-T08 vs T09-T10

df.cyc.age <- df5 %>% filter(time_point %in% c("T01", "T02", "T03", "T04", "T07", "T08", "T09", "T10")) %>%
  dplyr::select(-enrollment_status) %>% drop_na() %>%
  mutate(pre.post = ifelse(time_point %in% c("T01", "T02", "T07", "T08"), "pre", "post")) %>%
  mutate(storm = ifelse(time_point %in% c("T01", "T02", "T03", "T04"), "Cyclone Batsirai (February 2022)", "Cyclone Freddy (February 2023)")) %>%
  mutate(pre.post = factor(pre.post, levels = c("pre", "post"))) %>%
  group_by(storm, pre.post, age.cat.at.sample) %>%
  summarize(n = n(), n_pos = sum(rdt.result), prev = n_pos/n*100) %>% ungroup()  %>%
  mutate(age_cat = case_when(
    age.cat.at.sample == "A_0_2ys"      ~ "0 to 2 ys",
    age.cat.at.sample == "B_2_5ys"      ~ "2 to 5 ys",
    age.cat.at.sample == "C_6_13ys"     ~ "6 to 13 ys",
    age.cat.at.sample == "D_14_21ys"    ~ "14 to 21 ys",
    age.cat.at.sample == "E_22_40ys"    ~ "22 to 40 ys",
    age.cat.at.sample == "F_40plus_ys"  ~ "40 plus ys"
  ))

df.cyc.age <- df.cyc.age %>% mutate(age_cat = factor(age_cat, levels = c(unique(df.cyc.age$age_cat))))

df.cyc.age <- df.cyc.age %>% rowwise() %>% 
  mutate(Prev.CI.lower = f.CI_calc_lower(n_pos, n) * 100) %>%
  mutate(Prev.CI.upper = f.CI_calc_upper(n_pos, n) * 100)

p.cyc.age <- df.cyc.age %>% 
  ggplot(aes(x=pre.post, y=prev, fill=pre.post)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=Prev.CI.lower, ymax=Prev.CI.upper),
                linewidth = 0.25, # Thinner lines
                width = 0.15,       # Thinner horizontals on the Ts
                position="identity") +
  facet_grid(rows = vars(storm), cols = vars(age_cat)) +
  scale_fill_viridis_d(option = "mako", begin = 0.2, end = 0.9, name = "Cyclone\nperiod") +
  xlab("Comparing pre-cyclone and post-cyclone period") + 
  ylab("Percent of sampled individuals positive by RDT (%)") +
  labs(title    = "Malaria prevalence by age category",
       subtitle = "(all sites, 2 months before and after storm combined)") +
  theme_bw() + theme(legend.position = "right",
                     axis.title = element_text(size = 16),
                     strip.background = element_rect(fill = "white", colour = "white"),
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = "Andale mono"),
                     axis.text.y = element_text(size = 12),
                     strip.text = element_text(size = 12),
                     legend.text = element_text(size = 12), legend.title = element_text(size = 14))
p.cyc.age





