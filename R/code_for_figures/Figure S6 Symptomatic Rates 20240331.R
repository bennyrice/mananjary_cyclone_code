library(tidyverse)
library(patchwork)

#Reading in data (corresponding to function call in determining FOI script)
dfi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/Lab Projects/2020 Projects/CRS2020/1 DATA/DATA/DATA_QUESTIONNAIRE/R/output/cleaned_symptom_data 20240331.csv",
                       col_select = 1:24) %>%
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
  mutate(age.cat.at.sample.new = case_when(
    age.cat.at.sample == "A_0_2ys"     ~ "0-2ys",
    age.cat.at.sample == "B_2_5ys"     ~ "2-5ys",
    age.cat.at.sample == "C_6_13ys"    ~ "6-13ys",
    age.cat.at.sample == "D_14_21ys"   ~ "14-21ys",
    age.cat.at.sample == "E_22_40ys"   ~ "22-40ys",
    age.cat.at.sample == "F_40plus_ys" ~ "40+ys")) %>%
  mutate(age.cat.at.sample.new = factor(age.cat.at.sample.new, levels = c("0-2ys", "2-5ys", "6-13ys", "14-21ys", "22-40ys", "40+ys")))


#Calculating % with any fever within the last 2 weeks for RDT positives
df1 <- dfi %>% filter(rdt.result == 1) %>%
  mutate(fever.y.n = case_when(
    fever.recall.current.coded > 0 | fever.recall.2wks.coded > 0 ~ 1,
    .default = 0)) %>%
  group_by(age.cat.at.sample.new) %>% 
  summarize(n = n(),
            n.fever = sum(fever.y.n),
            perc.fever = sum(fever.y.n)/n()*100)

df1.overall <- dfi %>% filter(rdt.result == 1) %>%
  mutate(fever.y.n = case_when(
    fever.recall.current.coded > 0 | fever.recall.2wks.coded > 0 ~ 1,
    .default = 0)) %>%
  group_by(rdt.result) %>% 
  summarize(n = n(),
            n.fever = sum(fever.y.n),
            perc.fever = sum(fever.y.n)/n()*100)

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


#Adding confidence intervals
#Multiplying by 100 to make a percentage
df2 <- df1 %>% rowwise() %>% 
  mutate(CI.lower = f.CI_calc_lower(n.fever, n) * 100) %>%
  mutate(CI.upper = f.CI_calc_upper(n.fever, n) * 100)

df2.overall <- df1.overall %>% rowwise() %>% 
  mutate(CI.lower = f.CI_calc_lower(n.fever, n) * 100) %>%
  mutate(CI.upper = f.CI_calc_upper(n.fever, n) * 100)


#Plotting
p.fever <- df2 %>%
  ggplot(aes(x = age.cat.at.sample.new, y = perc.fever, fill = age.cat.at.sample.new)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x=age.cat.at.sample.new, ymin=CI.lower, ymax=CI.upper), 
                width=0.4, 
                alpha=0.9) +
  scale_fill_viridis_d(option = "viridis") +
  scale_y_continuous(limits = c(0, 65)) +
  xlab("Age category (years)") + ylab("Percent of RDT positive individuals reporting fever within the last 2 weeks") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", angle = 90, vjust = 0.5, hjust = 1))
p.fever

p.fever.overall <- df2.overall %>%
  ggplot(aes(x = rdt.result, y = perc.fever)) +
  geom_bar(stat = "identity", fill = "#AE123A") +
  geom_errorbar(aes(x=rdt.result, ymin=CI.lower, ymax=CI.upper), 
                width=0.4, 
                alpha=0.9) +
  scale_y_continuous(limits = c(0, 65)) +
  xlab("Overall\n(all ages)") + ylab("Percent of RDT positive individuals reporting fever within the last 2 weeks") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p.fever.overall
  
p.fever.overall + p.fever + plot_layout(widths = c(1, 6))


