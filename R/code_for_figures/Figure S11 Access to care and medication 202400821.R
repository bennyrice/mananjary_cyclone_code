library(tidyverse)
library(patchwork)

##############################################################################################################
## Reading in data
##############################################################################################################

#First, selecting columns to keep
cols.to.keep <- c(
  "full.code", "unique_ind_id", "site_code", "dob", "age_yrs", "age_cat", "sex", "enrollment_status", "time_point", 
  "rdt.result", 
  "sample.date", "age.yrs.at.sample", "age.cat.at.sample",
  "fever.recall.current", "fever.recall.current.coded", "fever.recall.2wks", "fever.recall.2wks.coded", 
  "rdt.recall", "rdt.recall.coded", "rdt.recall.other.loc", "rdt.recall.other.loc.coded", 
  
  "med.recall",
  "med.recall.anti.malaria.ACT",
  "med.recall.none",
  "med.recall.anti.bio.unsp",
  "med.recall.de.worm.unsp",
  "med.recall.iron.faf",
  "med.recall.anti.fever.unsp",
  "med.recall.anti.inflam.unsp",
  "med.recall.omepraz",
  "med.recall.anti.alg.unsp",
  "med.recall.anti.tb.tpi",
  "med.recall.vit.unsp",
  "med.recall.other",
  "med.recall.g.anti_malarials",
  "med.recall.g.none",
  "med.recall.g.anti_biotics",
  "med.recall.g.deworming",
  "med.recall.g.anti_anemia",
  "med.recall.g.pain_fever_reducer",
  "med.recall.g.gi_tract",
  "med.recall.g.anti_tb",
  "med.recall.g.vits_supps",
  "med.recall.g.other",
  
  "malaria.meds.y.n",
  "malaria.meds.y.n.coded",
  "malaria.meds",
  "malaria.meds.anti.malaria.ACT",
  "malaria.meds.anti.alg.para",
  "malaria.meds.pain_fever_reducer.uncoded",
  "malaria.meds.anti.bio.cotrim",
  "malaria.meds.anti.bio.amoxi",
  "malaria.meds.anti.malaria.quinine.cp",
  "malaria.meds.anti.malaria.quinine.inj",
  "malaria.meds.trad.aero.euc",
  "malaria.meds.trad.aero.orange",
  "malaria.meds.traditional.uncoded",
  "malaria.meds.western.other",
  "malaria.meds.trad.unsp",
  "malaria.meds.g.anti_malarials",
  "malaria.meds.g.pain_fever_reducer",
  "malaria.meds.g.anti_biotics",
  "malaria.meds.g.traditional",
  "malaria.meds.g.misc")


#Reading in data
dfi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/Lab Projects/2020 Projects/CRS2020/1 DATA/DATA/DATA_QUESTIONNAIRE/R/output/cleaned_symptom_data 20240331.csv",
                       col_select = all_of(cols.to.keep),
                       col_types = list(rdt.recall.other.loc = col_character(), rdt.recall.other.loc.coded = col_character())) %>%
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
  mutate(age.cat.at.sample.new = factor(age.cat.at.sample.new, levels = c("0-2ys", "2-5ys", "6-13ys", "14-21ys", "22-40ys", "40+ys"))) %>%
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
    time_point == "T11" ~ "T10"
  ))


##############################################################################################################
## Tidying the data
##############################################################################################################

df1 <- dfi %>% 
  #Coding RDT use
  mutate(rdt.recall.coded = case_when(
    rdt.recall.coded == "0"             ~ "No RDT performed",
    rdt.recall.coded == "1.csb.nearest" ~ "RDT at nearest CSB",
    rdt.recall.coded == "2.csb.other"   ~ "RDT at other CSB",
    rdt.recall.coded == "3.chw"         ~ "RDT by CHW",
    rdt.recall.coded == "4.hosp.mnj"    ~ "RDT at District Hospital",
    rdt.recall.coded == "5.priv.clin"   ~ "RDT by private provider",
  )) %>%
  mutate(rdt.recall.coded = factor(rdt.recall.coded, levels = rev(c("No RDT performed", "RDT at nearest CSB", "RDT at other CSB", "RDT by CHW", "RDT at District Hospital", "RDT by private provider")))) %>%
  #Coding meds: Anti-malarials, Anti-malarials & Others, Others, None
  mutate(mal.meds.dich = case_when(
    malaria.meds.g.anti_malarials == 1 &  malaria.meds.g.pain_fever_reducer == 0 & malaria.meds.g.anti_biotics == 0 & malaria.meds.g.traditional == 0 & malaria.meds.g.misc == 0  ~ "Anti-malarials",
    malaria.meds.g.anti_malarials == 1 & (malaria.meds.g.pain_fever_reducer == 1 | malaria.meds.g.anti_biotics == 1 | malaria.meds.g.traditional == 1 | malaria.meds.g.misc == 1) ~ "Anti-malarials and others",
    malaria.meds.g.anti_malarials == 0 & (malaria.meds.g.pain_fever_reducer == 1 | malaria.meds.g.anti_biotics == 1 | malaria.meds.g.traditional == 1 | malaria.meds.g.misc == 1) ~ "Others",
    is.na(malaria.meds.g.anti_malarials) ~ "None"
  )) %>%
  mutate(mal.meds.dich = factor(mal.meds.dich, levels = c("Anti-malarials", "Anti-malarials and others", "Others", "None")))


##############################################################################################################
## Diagnosis at baseline
##############################################################################################################
#Overall and by site: show proportion of RDT positive individuals who performed an RDT
df.dx.b <- df1 %>%   
  #Filtering to RDT positive individuals
  filter(rdt.result == 1) %>%
  #Filtering to baseline sample
  filter(time_point == "T0")

#number of samples
length(df.dx.b$rdt.result)

#relevant responses
table(df.dx.b$rdt.recall.coded)

p.dx.b <- df.dx.b %>% group_by(code.new, rdt.recall.coded) %>% summarize(n = n()) %>%
  ggplot(aes(x = code.new, y = n, fill = rdt.recall.coded)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis_d(name = "RDT use") +
  facet_wrap(vars(code.new), nrow = 1, scales = "free_x") +
  xlab("Site Code") +
  ylab("Number of infected individuals") +
  labs(title = "Use of malaria diagnosis at baseline (July-October 2021)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"))

##############################################################################################################
## Diagnosis thorought study
##############################################################################################################
#Overall and by site and time point: % of RDT positive individuals who performed an RDT outside the cohort
df.dx.c <- df1 %>% 
  #Filtering to RDT positive individuals
  filter(rdt.result == 1) %>%
  filter(time_point != "T0")

#number of samples
length(df.dx.c$rdt.result)

#relevant responses
table(df.dx.c$rdt.recall.coded)

p.dx.c <- df.dx.c %>% group_by(code.new, time_point, rdt.recall.coded) %>% summarize(n = n()) %>%
  ggplot(aes(x = time_point, y = n, fill = rdt.recall.coded)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis_d(name = "RDT use") +
  facet_wrap(vars(code.new), nrow = 1, scales = "free_x") +
  xlab("Time Point") +
  ylab("Number of infected individuals") +
  labs(title = "Use of malaria diagnosis during cohort study (2021-2023)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, family = "Andale Mono"),
        strip.background = element_rect(fill = "white"))


#Plotting both together
p.dx.b / p.dx.c

##############################################################################################################
## Treatment at baseline
##############################################################################################################
#Overall and by site: Porportion of RDT positive individuals who took meds
#Stacked bars by type of meds
df.rx.b <- df1 %>% 
  #Filtering to RDT positive individuals
  filter(rdt.result == 1) %>%
  filter(time_point == "T0")

#number of samples
length(df.rx.b$mal.meds.dich)

#relevant responses
table(df.rx.b$mal.meds.dich)

p.rx.b <- df.rx.b %>% group_by(code.new, mal.meds.dich) %>% summarize(n = n()) %>%
  ggplot(aes(x = code.new, y = n, fill = mal.meds.dich)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis_d(name = "Malaria treatment\nuse", option = "magma") +
  facet_wrap(vars(code.new), nrow = 1, scales = "free_x") +
  xlab("Site Code") +
  ylab("Number of infected individuals") +
  labs(title = "Use of malaria treatment at baseline (July-October 2021)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"))

##############################################################################################################
## Treatment thorought study
##############################################################################################################
#Overall and by site and time point: % of RDT positive individuals who took meds
#Stacked bars by type of meds
df.rx.c <- df1 %>% 
  #Filtering to RDT positive individuals
  filter(rdt.result == 1) %>%
  filter(time_point != "T0")

#number of samples
length(df.rx.c$mal.meds.dich)

#relevant responses
table(df.rx.c$mal.meds.dich)

p.rx.c <- df.rx.c %>% group_by(code.new, time_point, mal.meds.dich) %>% summarize(n = n()) %>%
  ggplot(aes(x = time_point, y = n, fill = mal.meds.dich)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis_d(name = "Malaria treatment\nuse", option = "magma") +
  facet_wrap(vars(code.new), nrow = 1, scales = "free_x") +
  xlab("Time Point") +
  ylab("Number of infected individuals") +
  labs(title = "Use of malaria treatment during cohort study (2021-2023)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, family = "Andale Mono"),
        strip.background = element_rect(fill = "white"))


#Plotting both together
p.rx.b / p.rx.c

(p.dx.b / p.dx.c)

(p.rx.b / p.rx.c)


##############################################################################################################
## Comparing baseline prevalence with access to diagnosis and treatment
##############################################################################################################

df2 <- df1 %>%
  filter(time_point == "T0") %>%
  #Calculating baseline prev
  group_by(code.new) %>% mutate(mal.prev = sum(rdt.result)/n()*100) %>% ungroup() %>%
  #Calculating baseline access to diagnosis
  group_by(code.new) %>% mutate(test.prev = sum(rdt.recall != "none")/n()*100) %>%
  #Calculating baseline access to treatment
  group_by(code.new, mal.prev, test.prev) %>% summarize(treat.prev = sum(malaria.meds.g.anti_malarials == 1, na.rm = TRUE)/n()*100)
  

p.test <- df2 %>%
  ggplot(aes(x = test.prev, y = mal.prev, label = code.new, color = code.new)) +
  geom_point() +
  scale_color_viridis_d(option = "turbo", name = "Site Code") +
  #xlim(0, 10) +
  xlab("Frequency of RDT testing\n(% of individuals reporting a prior RDT test)") +
  ylab("Prevalence of malaria at baseline (%)") +
  labs(title = "Comparing prevalence and diagnostic access") +
  theme_bw() +
  theme(panel.grid = element_blank())
p.test

p.treat <- df2 %>%
  ggplot(aes(x = treat.prev, y = mal.prev, label = code.new, color = code.new)) +
  geom_point() +
  scale_color_viridis_d(option = "turbo", name = "Site Code") +
  #xlim(0, 10) +
  xlab("Frequency of treatment\n(% of individuals reporting prior use of anti-malarials)") +
  ylab("Prevalence of malaria at baseline (%)") +
  labs(title = "Comparing prevalence and treatment access") +
  theme_bw() +
  theme(panel.grid = element_blank())
p.treat

p.test + p.treat

