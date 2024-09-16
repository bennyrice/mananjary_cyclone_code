library(tidyverse)
library(patchwork)

#Reading in data (corresponding to function call in determining FOI script)
dfi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/Lab Projects/2020 Projects/CRS2020/1 DATA/DATA/DATA_BEDNET/bednet_data_20230907.csv")

df.rdt <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_malaria_rdt/rdt_data_20231007.csv")

df1 <- full_join(df.rdt, dfi, by = join_by(unique_ind_id, time_point, site_code)) %>%
  filter(!is.na(rdt.result)) %>% filter(!is.na(bednet.24hrs.coded)) %>%
  mutate(code.new = case_when(
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
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10"))) %>%
  mutate(time_point = case_when(
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
    time_point == "T11" ~ "T10")) %>%
  #24hr recall
  mutate(bednet.stat.24hr = case_when(
    bednet.condition.24hrs == "good"        ~ "Yes",
    bednet.condition.24hrs == "bad_other"   ~ "Yes (damaged)",
    bednet.condition.24hrs == "hole_finger" ~ "Yes (damaged)",
    bednet.condition.24hrs == "hole_fist"   ~ "Yes (damaged)",
    bednet.condition.24hrs == "hole_head"   ~ "Yes (damaged)",
    bednet.condition.24hrs == "hole_unsp"   ~ "Yes (damaged)",
    is.na(bednet.condition.24hrs)           ~ "No")) %>%
  mutate(bednet.stat.24hr = factor(bednet.stat.24hr, levels = c("No", "Yes (damaged)", "Yes"))) %>%
  #2wk recall
  mutate(bednet.stat.2wks = case_when(
    bednet.condition.2wks == "good"        ~ "Yes",
    bednet.condition.2wks == "bad_other"   ~ "Yes (damaged)",
    bednet.condition.2wks == "hole_finger" ~ "Yes (damaged)",
    bednet.condition.2wks == "hole_fist"   ~ "Yes (damaged)",
    bednet.condition.2wks == "hole_head"   ~ "Yes (damaged)",
    bednet.condition.2wks == "hole_unsp"   ~ "Yes (damaged)",
    is.na(bednet.condition.2wks)           ~ "No")) %>%
  mutate(bednet.stat.2wks = factor(bednet.stat.2wks, levels = c("No", "Yes (damaged)", "Yes")))
  



#Due to sampling schedule, looking by month is not helpful as biased by some sites not covered during certain months
df2 <- df1 %>% mutate(month = lubridate::floor_date(sample.date, "month")) %>%
  group_by(month, bednet.stat.24hr) %>% 
  summarize(n = n()) %>%
  mutate(prop = n / sum(n)) %>% ungroup()

p.bednet.month <- df2 %>%
  ggplot(aes(x = month, y = prop, fill = bednet.stat.24hr)) +
  geom_bar(position = "fill", stat = "identity")
p.bednet.month


#by time point (all sites combined)
df3 <- df1 %>% 
  group_by(time_point, bednet.stat.24hr) %>% 
  summarize(n = n()) %>%
  mutate(prop = n / sum(n)) %>% ungroup()

p.bednet.tp <- df3 %>%
  ggplot(aes(x = time_point, y = prop, fill = bednet.stat.24hr)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Bednet Usage") +
  xlab("Sampling Time Point") + ylab("Proportion of individuals reporting using a bednet in last 24 hours") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5))
p.bednet.tp







#by site and time point
df4 <- df1 %>% 
  dplyr::select(all_of(c("code.new", "unique_ind_id", "sample.date", "time_point", "bednet.stat.24hr", "bednet.stat.2wks"))) %>%
  pivot_longer(!code.new:time_point, names_to = "period", values_to = "status") %>%
  group_by(code.new, time_point, period, status) %>% 
  summarize(n = n()) %>%
  mutate(prop = n / sum(n)) %>% ungroup()

p.bednet.tp.site <- df4 %>%
  ggplot(aes(x = time_point, y = prop, fill = status)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(cols = vars(code.new), rows = vars(period)) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Bednet Usage") +
  xlab("Sampling Time Point") + ylab("Proportion of individuals reporting using a bednet") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5),
        strip.background = element_rect(fill = "white"))
p.bednet.tp.site

#2wk and 24hr recall are similar
df4 %>%
  ggplot(aes(x = period, y = prop, fill = status)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(vars(time_point), nrow = 1) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Bednet Usage") +
  xlab("period") + ylab("Proportion of individuals reporting using a bednet") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5),
        strip.background = element_rect(fill = "white"))

##Only showing 24hr since similar to 2wks
##Marking cyclone period (February 2022 and February 2023)
p.bednet.tp.site.24hr <- df4 %>%
  filter(period == "bednet.stat.24hr") %>%
  ggplot(aes(x = time_point, y = prop, fill = status)) +
  geom_bar(position = "fill", stat = "identity") +
  #Adding a vertical line for storm impact
  #geom_vline(xintercept = c(3, 10), linetype = 2, alpha = 0.9, color = "grey20") +
  facet_wrap(vars(code.new), nrow = 1) +
  scale_fill_viridis_d(option = "cividis", direction = -1, name = "Bednet Usage\n(24hr Recall)") +
  xlab("Sampling Time Point") + ylab("Proportion of individuals reporting using a bednet\n(in previous 24 hours)") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5),
        strip.background = element_rect(fill = "white"))
p.bednet.tp.site.24hr

##Prevalence heat strip: By time point
p.prev.heat <- df1 %>%
  group_by(code.new, time_point) %>% summarize(prev = sum(rdt.result)/n()*100) %>%
  ggplot(aes(x = time_point, y = 1, fill = prev)) +
  geom_tile() +
  facet_wrap(vars(code.new), nrow = 1) +
  scale_fill_viridis_c(option = "inferno", name = "Prevalence\n(%)", limits = c(0, 65)) +
  ylab(NULL) +
  xlab("Prevalence by time point") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = "white"))
p.prev.heat

p.prev.heat.avg <- df1 %>%
  group_by(code.new) %>% summarize(prev = sum(rdt.result)/n()*100) %>%
  ggplot(aes(x = code.new, y = 1, fill = prev)) +
  geom_tile() +
  facet_wrap(vars(code.new), nrow = 1, scales = "free_x") +
  scale_fill_viridis_c(option = "inferno", name = "Avg. Prevalence\n(%)", limits = c(0, 65)) +
  ylab(NULL) +
  xlab("Average prevalence") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = "white"))
p.prev.heat.avg

#Adding sample size
p.n <- df1 %>%
  group_by(code.new, time_point) %>% summarize(n = n()) %>%
  ggplot(aes(x = time_point, y = n)) +
  geom_point() +
  facet_wrap(vars(code.new), nrow = 1) +
  ylim(0, 375) +
  ylab("n") +
  xlab("Sample size by time point (number of individuals sampled)") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5),
        strip.background = element_rect(fill = "white"))
p.n

p.bednet.tp.site.24hr / p.n / p.prev.heat / p.prev.heat.avg + plot_layout(heights = c(10, 2, 1, 1))


################################################################################################
#pre-post storms with confidence intervals

#Checking dates
df1 %>% group_by(code.new, time_point, sample.date) %>% summarize(n = n()) %>%
  ggplot(aes(x = sample.date, y = n, color = code.new)) +
  geom_point() +
  scale_color_viridis_d(option = "turbo") +
  facet_grid(cols = vars(time_point), rows = vars(code.new), scales = "free")

#Batsirai
#Sites:
# "MNJ.01" Pre: T01 ; Post: T02
# "MNJ.02" Pre: T01 ; Post: T02
# "MNJ.03" Pre: T01 ; Post: T02
# "MNJ.04" Pre: T01 ; Post: T02
# "MNJ.05" Pre: T01 ; Post: T03
# "MNJ.06" Pre: T02 ; Post: T03
# "MNJ.07" Pre: T02 ; Post: T03
# "MNJ.08" Pre: T02 ; Post: T03
# "MNJ.09" Pre: T02 ; Post: T03
# "MNJ.10" Pre: T02 ; Post: T03

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

#Cyclone Batsirai: T02-T01 vs T04-T03
#Cyclone Freddy:   T09-T08 vs T11-T10
p.T01.T03 <- df1 %>% filter(time_point %in% c("T01", "T03")) %>%
  group_by(time_point, bednet.stat.24hr) %>% 
  summarize(n = n()) %>%
  mutate(prop = n / sum(n)) %>% ungroup() %>%
  #Calculating the proportion not using a bednet
  group_by(time_point) %>% mutate(sum_n = sum(n)) %>% filter(bednet.stat.24hr == "No") %>%
  #Adding CIs
  rowwise() %>% 
  mutate(CI.lower = f.CI_calc_lower(n, sum_n)) %>%
  mutate(CI.upper = f.CI_calc_upper(n, sum_n)) %>%
  #Plotting
  ggplot(aes(x = time_point, y = prop)) +
  geom_bar(stat = "identity", fill = "#FCEA66") +
  geom_errorbar(aes(x=time_point, ymin=CI.lower, ymax=CI.upper), 
                width=0.4, 
                alpha=0.9) +
  #Adding a vertical line for storm impact
  geom_vline(xintercept = 1.5, linetype = 2, alpha = 0.9, color = "grey20") +
  scale_y_continuous(limits = c(0, 0.5)) +
  xlab("Sampling Time Point") + ylab("Proportion of individuals not using a bednet\n(in previous 24 hours)") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5))
p.T01.T03


p.T08.T10 <- df1 %>% filter(time_point %in% c("T08", "T10")) %>%
  group_by(time_point, bednet.stat.24hr) %>% 
  summarize(n = n()) %>%
  mutate(prop = n / sum(n)) %>% ungroup() %>%
  group_by(time_point) %>% mutate(sum_n = sum(n)) %>% filter(bednet.stat.24hr == "No") %>%
  #Adding CIs
  rowwise() %>% 
  mutate(CI.lower = f.CI_calc_lower(n, sum_n)) %>%
  mutate(CI.upper = f.CI_calc_upper(n, sum_n)) %>%
  ggplot(aes(x = time_point, y = prop)) +
  geom_bar(stat = "identity", fill = "#FCEA66") +
  geom_errorbar(aes(x=time_point, ymin=CI.lower, ymax=CI.upper), 
                width=0.4, 
                alpha=0.9) +
  #Adding a vertical line for storm impact
  geom_vline(xintercept = 1.5, linetype = 2, alpha = 0.9, color = "grey20") +
  scale_y_continuous(limits = c(0, 0.5)) +
  xlab("Sampling Time Point") + ylab("Proportion of individuals not using a bednet\n(in previous 24 hours)") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0.5))
p.T08.T10

p.T01.T03 + p.T08.T10





