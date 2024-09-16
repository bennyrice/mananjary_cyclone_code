library(tidyverse)
library(patchwork)



## Reading in data ##################################################
df <- read.csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/mnj_rdt_data_foi_20230821.csv")

## Making individual ID code and household ID code a factor for fitting random effects
df$unique_ind_id <- as.factor(df$unique_ind_id)
## First 5 digits of the individual ID code give the site code + household code
df$household <- as.factor(substring(df$unique_ind_id,1,5)) 

################################################################################
## Sensitivty analyses
################################################################################

#For each loop:
#   Simulate RDT results using rbinom
#   Fit the gam
#   Predict Pr(inf) after 1 day for each site and a standard age (10 yo)
#   Store in a data frame for plotting

f.sensitivity <- function(data, sensitivity, specificity, site, reps){
  
  ## Step 1: Simulated data: Flipping RDT results based on 90% accuracy
  
  # 90% sensitivity = approx. 10% of negatives are false negatives --> flip 10% of negatives to positives
  # 90% specificity = approx. 10% of positives are false positives --> flip 10% of positives to negatives
  
  # If test negative, draw a 1 with 10% probability and a 0 with 90% probability etc
  
  #Creating empty vectors to hold data
  v.trial        <- rep(NA, reps)
  v.n_pos        <- rep(NA, reps)
  prob.inf.trial <- rep(NA, reps)
  #Splitting the data frame to quicken code
  df.sa <- data %>% mutate(rdt.result.sa = rdt.result) %>% 
    select(unique_ind_id, site_code, age_yrs, interlude, interlude.days, rdt.result)
  df.sa.0 <- df.sa %>% filter(rdt.result == 0)
  df.sa.1 <- df.sa %>% filter(rdt.result == 1)
  
  for(i in 1:reps){
    
    #Splitting then sampling/simulating then re-binding
    df.sa.0.i <- df.sa.0 %>% mutate(rdt.result.sa = rbinom(length(rdt.result), 1, 1-sensitivity))
    df.sa.1.i <- df.sa.1 %>% mutate(rdt.result.sa = rbinom(length(rdt.result), 1, specificity))
    df.sa.i <- rbind(df.sa.0.i, df.sa.1.i) %>% arrange(unique_ind_id, interlude)
    
    v.trial[i] <- i
    v.n_pos[i] <- sum(df.sa.i$rdt.result.sa)
    
    fit.sa <- gam(rdt.result.sa~s(age_yrs)+site_code+offset(log(interlude.days)),
                  family=binomial(link="cloglog"), data=df.sa.i)
    
    newData.sa <- data.frame(age_yrs        = 10,
                             site_code      = site,
                             interlude.days = 1)
    
    tmp.sa <- predict(fit.sa, newdata=newData.sa, se=FALSE, type="response")
    prob.inf.trial[i] <- tmp.sa
    
  }
  
  df.sa <- tibble(v.trial, v.n_pos, prob.inf.trial)
  return(df.sa)
}


df.sa1 <- f.sensitivity(data = df, sensitivity = 1.0, specificity = 1.0, site = "S3", reps = 1)
df.sa2 <- f.sensitivity(data = df, sensitivity = 0.9, specificity = 0.9, site = "S3", reps = 1)
df.sa3 <- f.sensitivity(data = df, sensitivity = 0.8, specificity = 0.8, site = "S3", reps = 1)


start.time <- Sys.time()
#
df.sa1 <- f.sensitivity(data = df, sensitivity = 0.9, specificity = 0.9, site = "S3", reps = 50)
#
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken


v.sites <- c("S1", "S2", "S3", "S5", "S6", "N1", "N2", "N3", "N4", "N5")

sa.list.100.100 <- list(NA)
sa.list.99.99   <- list(NA)
sa.list.95.95   <- list(NA)
sa.list.90.90   <- list(NA)

for(i in 1:length(v.sites)){
  sa.list.100.100[[i]]   <- f.sensitivity(data = df, sensitivity = 1.00, specificity = 1.00, site = v.sites[i], reps = 1) %>%
    mutate(site_code = v.sites[i])
}
for(i in 1:length(v.sites)){
  sa.list.99.99[[i]]     <- f.sensitivity(data = df, sensitivity = 0.99, specificity = 0.99, site = v.sites[i], reps = 500) %>%
    mutate(site_code = v.sites[i])
}
for(i in 1:length(v.sites)){
  sa.list.95.95[[i]]     <- f.sensitivity(data = df, sensitivity = 0.95, specificity = 0.95, site = v.sites[i], reps = 500) %>%
    mutate(site_code = v.sites[i])
}
for(i in 1:length(v.sites)){
  sa.list.90.90[[i]]     <- f.sensitivity(data = df, sensitivity = 0.90, specificity = 0.90, site = v.sites[i], reps = 500) %>%
    mutate(site_code = v.sites[i])
}


sa.100.100 <- bind_rows(sa.list.100.100) %>% mutate(sens.spec = "100%, 100%")
sa.99.99   <- bind_rows(sa.list.99.99)   %>% mutate(sens.spec = "99%, 99%")
sa.95.95   <- bind_rows(sa.list.95.95)   %>% mutate(sens.spec = "95%, 95%")
sa.90.90   <- bind_rows(sa.list.90.90)   %>% mutate(sens.spec = "90%, 90%")



df.sensitivity <- rbind(sa.100.100, sa.99.99, sa.95.95, sa.90.90) %>%
  mutate(sens.spec = factor(sens.spec, levels = c("100%, 100%", 
                                                  "99%, 99%", 
                                                  "95%, 95%", 
                                                  "90%, 90%"))) %>%
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
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10")))


#export to avoid refitting
# write_csv(df.sensitivity, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/sensitivity_analyses_90_90.csv")

df.sensitivity.means <- df.sensitivity %>% group_by(code.new, sens.spec) %>% summarize(mean = mean(prob.inf.trial))

p.sensitivity <- df.sensitivity %>%
  ggplot(aes(x = prob.inf.trial, fill = sens.spec, color = sens.spec)) +
  #*change binwidth
  geom_histogram(binwidth = 0.00005, position = "identity", alpha = 0.5, linewidth = 0.2) +
  #Adding vertical lines for the means for the scenarios
  geom_vline(data = df.sensitivity.means,
             aes(xintercept = mean, color = sens.spec),
             linetype = 2, alpha = 0.7) +
  #Facet wrapping: 2 x 5
  facet_wrap(vars(code.new), nrow = 2) +
  #*check for faint yellow at end of scale
  scale_fill_viridis_d(option = "viridis", name = "Sensitivity and\nSpecificity\nScenario", end = 1) +
  scale_color_viridis_d(option = "viridis", name = "Sensitivity and\nSpecificity\nScenario", end = 1) +
  xlab("Estimated Force of Infection (FOI)") + ylab("Count (from 500 simulations)") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono"),
        strip.background = element_rect(fill = "white"))
p.sensitivity







