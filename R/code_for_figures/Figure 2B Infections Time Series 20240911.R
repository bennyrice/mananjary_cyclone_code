library(tidyverse)
library(patchwork)

#read in RDT data
dfi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_malaria_rdt/rdt_data_20231007.csv") 

#### Summary Stats: RDT sampling
df.mean.sampling <- dfi %>% 
  #drop missing data
  filter(!is.na(rdt.result)) %>%  
  #count number of observations per individual
  group_by(unique_ind_id) %>% summarize(n = n())
# Median and mean number of samples per individual
mean(df.mean.sampling$n)
median(df.mean.sampling$n)
# Number of individuals for a given sample completeness
length(df.mean.sampling$n[df.mean.sampling$n > 0])
length(df.mean.sampling$n[df.mean.sampling$n > 1])
length(df.mean.sampling$n[df.mean.sampling$n > 2])
length(df.mean.sampling$n[df.mean.sampling$n > 3])
length(df.mean.sampling$n[df.mean.sampling$n > 4])
length(df.mean.sampling$n[df.mean.sampling$n > 5])
length(df.mean.sampling$n[df.mean.sampling$n > 6])
length(df.mean.sampling$n[df.mean.sampling$n > 7])

#### Summary Stats: Time interval
#read in RDT data
df.int <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/mnj_rdt_data_foi_20230821.csv") 
# Median and mean length of interval between samples
mean(df.int$interlude.days)
median(df.int$interlude.days)
hist(df.int$interlude.days)

#clean up RDT data for plotting
df1 <- dfi %>%
  #drop rows without data
  filter(!is.na(rdt.result)) %>%
  #drop baseline time_point
  filter(time_point != "T01") %>%
  #filter out individuals with >50% missing data
  group_by(unique_ind_id) %>% mutate(n = n()) %>% ungroup() %>% filter(n > 6)
#randomly sample 100 individuals from each site
df.sample <- df1 %>% group_by(site_code, unique_ind_id) %>% summarize(n = n()) %>% ungroup() %>%
  group_by(site_code) %>% slice_sample(n = 100)
#subset to the randomly sampled 100
df2 <- df1 %>% filter(unique_ind_id %in% df.sample$unique_ind_id) 
#cleaning up labels for plotting
df3 <- df2 %>%
  #moving baseling time point to t0 and subsequent follow up time points to T01-T10
  mutate(time_point = case_when(
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
  #providing standardized site codes
  mutate(code.new = case_when(
    site_code == "N5" ~ "MNJ.01",
    site_code == "N4" ~ "MNJ.02",
    site_code == "N3" ~ "MNJ.03",
    site_code == "N2" ~ "MNJ.04",
    site_code == "N1" ~ "MNJ.05",
    site_code == "S1" ~ "MNJ.06",
    site_code == "S2" ~ "MNJ.07",
    site_code == "S3" ~ "MNJ.08",
    site_code == "S5" ~ "MNJ.09",
    site_code == "S6" ~ "MNJ.10")) %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06",
                                                "MNJ.04",
                                                "MNJ.03",
                                                "MNJ.05",
                                                "MNJ.02",
                                                "MNJ.07",
                                                "MNJ.08",
                                                "MNJ.09",
                                                "MNJ.01",
                                                "MNJ.10")))


#Plotting
p.2b <- df3 %>% ggplot(aes(x = time_point, y = unique_ind_id, color = factor(rdt.result))) +
  geom_point(size = 1.5) +
  facet_wrap(vars(code.new), scales = "free_y", nrow = 1) +
  scale_color_manual(name = "Malaria status\n(by RDT)",
                     values = c("grey20", "red"),
                     labels = c("Negative", "Positive")) +
  scale_x_discrete(breaks = c("T01", "T10")) +
  xlab("Time Point") + ylab("Individual") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 0),
        axis.title = element_text(family = "Arial", size = 16),
        panel.grid.major.y = element_line(linewidth = 0.2, color = "grey60"),
        panel.grid.major.x = element_blank())
p.2b
