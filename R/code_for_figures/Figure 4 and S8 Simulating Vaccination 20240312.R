library(tidyverse)
library(ggrepel)

####################################################################################################
# Defining probabilities
####################################################################################################

#Defining parameters and using dummy values to valide equations
#Probability of full vaccination (v)
v   <- 0.50
#Probability of partial vaccination (v_p)
v_p   <- 0
#Probability of infection (p_I)
p_I <- 0.45
#Vaccine efficacy against infection (E_i)
E_i <- 0.04
#Vaccine efficacy against symptomatic infection (E_s)
E_s <- 0.80
#Average efficacy for partial vaccinated individuals (E_p)
E_p <- 0.40
#Probability of clinical symptoms (ie fever) upon infection (s)
s   <- 0.53


#Model 1: No partial vaccination

#1 Symptomatic infections among vaccinated individuals (I_V_S)
I_V_S = v * p_I * (1 - E_i) * s * (1 - E_s)
#2 Asymptomatic infections among vaccinated individuals (I_V_nS)
I_V_nS = v * p_I * (1 - E_i) * (1 - s * (1 - E_s))
#3 Uninfected vaccinated individuals (nI_V)
nI_V = v * (1 - (p_I * (1 - E_i)))
#4 Symptomatic infections among unvaccinated individuals (I_nV_S)
I_nV_S = (1 - v) * p_I * s
#5 Asymptomatic infections among unvaccinated individuals (I_nV_nS)
I_nV_nS = (1 - v) * p_I * (1 - s)
#6 Uninfected unvaccinated individuals (nI_nV)
nI_nV = (1 - v) * (1 - p_I)

#Verifying that sum of probabilities for outcomes is 1
I_V_S + I_V_nS + nI_V + I_nV_S + I_nV_nS + nI_nV


#Model 1: No partial vaccination

#1 Symptomatic infections among fully vaccinated individuals (I_V_S)
I_V_S = v * p_I * (1 - E_i) * s * (1 - E_s)
#2 Asymptomatic infections among fully vaccinated individuals (I_V_nS)
I_V_nS = v * p_I * (1 - E_i) * (1 - s * (1 - E_s))
#3 Uninfected fully vaccinated individuals (nI_V)
nI_V = v * (1 - (p_I * (1 - E_i)))

#4 Symptomatic infections among partially vaccinated individuals (I_V_S)
I_VP_S = v_p * p_I * (1 - E_i) * s * (1 - E_p)
I_VP_S
#5 Asymptomatic infections among partially vaccinated individuals (I_V_nS)
I_VP_nS = v_p * p_I * (1 - E_i) * (1 - s * (1 - E_p))
I_VP_nS
#6 Uninfected fully vaccinated individuals (nI_V)
nI_VP = v_p * (1 - (p_I * (1 - E_i)))

#7 Symptomatic infections among unvaccinated individuals (I_nV_S)
I_nV_S = (1 - v - v_p) * p_I * s
#8 Asymptomatic infections among unvaccinated individuals (I_nV_nS)
I_nV_nS = (1 - v - v_p) * p_I * (1 - s)
#9 Uninfected unvaccinated individuals (nI_nV)
nI_nV = (1 - v - v_p) * (1 - p_I)

#Verifying that sum of probabilities for outcomes is 1
I_V_S + I_V_nS + nI_V + I_VP_S + I_VP_nS + nI_VP + I_nV_S + I_nV_nS + nI_nV



####################################################################################################
# Calculating symptomatic infection rates under vaccination scenarios
####################################################################################################

# Symptomatic infections = fully vaccinated symptomatics (I_V_S) + partially vaccinated symptomatics (I_VP_S) + unvaccinated symptomatics (I_nV_S)
# I_V_S + I_VP_S + I_nV_S = v * p_I * (1 - E_i) * s * (1 - E_s) + v_p * p_I * (1 - E_i) * s * (1 - E_p) + (1 - v - v_p) * p_I * s

### Test case

#Params
v <- 0.7
v_p <- 0
p_I <- 0.25
E_i <- 0
E_s <- 0.55
E_p <- 0
s <- 0.4


#Probability for test parameter set
PR.sympt.inf <- v * p_I * (1 - E_i) * s * (1 - E_s) + v_p * p_I * (1 - E_i) * s * (1 - E_p) + (1 - v - v_p) * p_I * s
PR.sympt.inf

#First, building a data frame for a range of vaccination coverages, symptomatic rates, and vaccine efficacies
#Fixing probability of infection to 10% (since looking at proportions only)
#No partial vaccination or efficacy against infection yet

#paramater ranges:
# p_I = 0.1
# v: 0:100
# s: 1:100
# E_s: 0:100
#n = approx. 101^3 = approx. 1,030,301 combinations for a given level of infection

#Making an initial data frame for vaccination simulation
#using expand.grid() to write out all possible combinations
dfv.i <- tibble(expand.grid(v = (0:100)/100,
                            s = (1:100)/100, #no cases expected if symptomatic infection = 0
                            E_s = (0:100)/100))

dfv <- dfv.i %>%
  #fixing the probability of infection 0.1 here
  mutate(p_I = 0.1) %>%
  #setting fixed paramaters
  mutate(E_i = 0, v_p = 0, E_p = 0) %>%
  #Calculate the probability of a symptomatic infection
  mutate(PR.sympt.inf = v * p_I * (1 - E_i) * s * (1 - E_s) + v_p * p_I * (1 - E_i) * s * (1 - E_p) + (1 - v - v_p) * p_I * s) %>%
  #Calculating baseline number of expected cases and declines under vaccination
  mutate(baseline_proportion = p_I*s,
         perc_decline = (baseline_proportion - PR.sympt.inf)/baseline_proportion*100,
         cases_averted_p1000 = baseline_proportion*1000 - PR.sympt.inf*1000,
         baseline_case_count_p1000 = baseline_proportion*1000,
         vaccinated_case_count_p1000 = PR.sympt.inf*1000) %>%
  #showing vaccinated proportion as a percentage for plotting
  mutate(v_perc = v*100)

##################################################################
## Plotting % reduction in cases
##################################################################

p1 <- dfv %>% 
  #Fixing symptomatic rate at an arbitrary level to reduce data frame size (s has no effect here)
  filter(s == 0.28) %>% 
  ggplot(aes(x = v_perc, y = perc_decline, group = E_s, color = E_s)) +
  #Shading portion of the plot to higlight v and E combos with reduction > 50%
  annotate("rect", xmin=50, xmax=100, ymin=50, ymax=100, fill = "black", alpha = 0.07) +
  #Adding the geom_path
  geom_path(linewidth = 0.8, alpha = 0.8) +
  #Overlay R21
  geom_path(data = dfv %>% filter(s == 0.28) %>% filter(E_s == 0.75), 
            aes(x = v_perc, y = perc_decline, group = E_s),
            color = "black", linewidth = 1.3, alpha = 0.3) +
  geom_label_repel(data = dfv %>% filter(s == 0.28) %>% filter(E_s == 0.75) %>% filter(v == 1), 
                   aes(label = paste0("E_s = ", E_s, " (R21)")),
                   nudge_x = 6.1,
                   na.rm = TRUE) +
  #Overlay RTS,S
  geom_path(data = dfv %>% filter(s == 0.28) %>% filter(E_s == 0.55), 
            aes(x = v_perc, y = perc_decline, group = E_s),
            color = "black", linewidth = 1.3, alpha = 0.3) +
  geom_label_repel(data = dfv %>% filter(s == 0.28) %>% filter(E_s == 0.55) %>% filter(v == 1), 
                   aes(label = paste0("E_s = ", E_s, " (RTS,S)")),
                   nudge_x = 7.2,
                   na.rm = TRUE) +
  #Adding horizontal and vertical guidelines
  #50% reduction threshold
  geom_segment(aes(x=0, xend=100, y=50, yend=50), color="grey50", linetype="dashed", linewidth=0.3) +
  #R21
  geom_segment(aes(x=50/0.75, xend=50/0.75, y=50, yend=100), color="grey50", linetype="dashed", linewidth=0.3) +
  #RTSS
  geom_segment(aes(x=50/0.55, xend=50/0.55, y=50, yend=100), color="grey50", linetype="dashed", linewidth=0.3) +
  #Perfect vaccine
  geom_segment(aes(x=50/1.00, xend=50/1.00, y=50, yend=100), color="grey50", linetype="dashed", linewidth=0.3) +
  #Overlaying intersection points
  #R21
  geom_point(aes(x = 50/0.75, y = 50), color = "black", alpha = 0.6) +
  #RTS,S
  geom_point(aes(x = 50/0.55, y = 50), color = "black", alpha = 0.6) +
  #Perfect vaccine
  geom_point(aes(x = 50/1.00, y = 50), color = "black", alpha = 0.6) +
  #Scales and themes
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_reverse(breaks = seq(0, 100, by = 10)) +
  scale_color_viridis_c(option = "cividis",
                        name = "Vaccine\nEfficacy",
                        end = 1) +
  xlab(expression(paste("Percent fully vaccinated at time of disruptive event (", italic("v"), ")"))) +
  ylab(expression(paste("Percent reduction in symptomatic infections expected (", italic("vE"), ")"))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16, family = "Arial"),
        axis.text = element_text(size = 13, family = "Arial"))
p1

#Pulling out some bounds
#Symptomatic rates for children: 37.0-53.2%
#R21
#min scenario
dfv %>% filter(s == 0.370) %>% filter(E == 0.68) %>% filter(v == 0.50)
dfv %>% filter(s == 0.370) %>% filter(E == 0.75) %>% filter(v == 0.70)
#max scenario
dfx <- dfv %>% filter(s == 0.530) %>% filter(E == 0.68) %>% filter(v == 0.70)
dfx <- dfv %>% filter(s == 0.530) %>% filter(E == 0.75) %>% filter(v == 0.70)

#RTSS
dfv %>% filter(s == 0.28) %>% filter(E == 0.50) %>% filter(v == 0.50)
dfv %>% filter(s == 0.28) %>% filter(E == 0.56) %>% filter(v == 0.70)


##################################################################
## Plotting hazard ratio time series
##################################################################

#Fix coverage at 70%, s at 50%
#FOI: age group: 5 yos
#High FOI site: MNJ.10
#Mid  FOI site: MNJ.07
#Low  FOI site: MNJ.06
#Start date: Freddy landfall date (February 21, 2023)

#FOI time series from foi_post_storm but run out to 90 days (ie 3 months)


df.inf <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_post_storm.csv") %>%
  filter(storm == "CYCLONE FREDDY", age == 5) %>%
  #Subsetting to high - mid - low exemplar sites for plotting
  filter(code.new %in% c("MNJ.10", "MNJ.07", "MNJ.06")) %>%
  mutate(FOI_cat = case_when(
    code.new == "MNJ.10" ~ "foi.high",
    code.new == "MNJ.07" ~ "foi.mid",
    code.new == "MNJ.06" ~ "foi.low"))

dfv.hazard <- df.inf %>% 
  dplyr::select(day, prob.inf, FOI_cat) %>%
  pivot_wider(values_from = prob.inf, names_from = FOI_cat) %>%
  mutate(v = 0.70,
         s = 0.50,
         #E.90 = 0.90,
         #E.80 = 0.80,
         E.75 = 0.75, #R21
         #E.70 = 0.70,
         #E.60 = 0.60,
         E.55 = 0.55, #RTS,S
         #E.50 = 0.50,
         #E.40 = 0.40,
         #E.30 = 0.30,
         #E.20 = 0.20,
         #E.10 = 0.10,
         E.0 = 0) %>%
  pivot_longer(!c(day:s), values_to = "E_s", names_to = "efficacy.cat") %>%
  pivot_longer(!c(day, v, s, efficacy.cat, E_s), values_to = "p_I", names_to = "foi.cat") %>%
  #setting fixed paramaters
  mutate(E_i = 0, v_p = 0, E_p = 0) %>%
  #Calculate the probability of a symptomatic infection
  mutate(PR.sympt.inf = v * p_I * (1 - E_i) * s * (1 - E_s) + v_p * p_I * (1 - E_i) * s * (1 - E_p) + (1 - v - v_p) * p_I * s) %>%
  #Calculating baseline number of expected cases and declines under vaccination
  mutate(baseline_proportion = p_I*s,
         perc_decline = (baseline_proportion - PR.sympt.inf)/baseline_proportion*100,
         cases_averted_p1000 = baseline_proportion*1000 - PR.sympt.inf*1000,
         baseline_case_count_p1000 = baseline_proportion*1000,
         vaccinated_case_count_p1000 = PR.sympt.inf*1000) %>%
  #showing vaccinated proportion as a percentage for plotting
  mutate(v_perc = v*100) %>%
  #cleaning up factors
  mutate(foi.cat = factor(foi.cat, levels = rev(c("foi.high", "foi.mid", "foi.low"))))

foi.cat_names <- c(
  foi.high = "High FOI locality (MNJ.10)",
  foi.mid  = "Mid FOI locality (MNJ.07)",
  foi.low  = "Low FOI locality (MNJ.06)"
)

p2 <- dfv.hazard %>% ggplot(aes(x = day, y = PR.sympt.inf, group = efficacy.cat)) +
  geom_line(aes(linetype = efficacy.cat), linewidth = 1, color = "#AE123A") +
  facet_wrap(vars(foi.cat), labeller = labeller(foi.cat = foi.cat_names)) +
  scale_linetype_discrete(name = "Vaccination Scenario",
                          labels = c("No vaccine", "Es = 55%", "Es = 75%")) +
  xlab("Days after storm landfall") + 
  ylab("Probability of symptomatic malaria infection (among young children)") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        #legend.position = "bottom",
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12))
p2


##################################################################
## Plotting cases averted
##################################################################

dfc.i <- tibble(expand.grid(v = (0:100)/100,
                            E_s = (0:100)/100,
                            p_I = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                            #no cases expected if symptomatic infection = 0
                            s = c(0.3, 0.4, 0.5, 0.6)))

dfc <- dfc.i %>%
  #setting fixed paramaters
  mutate(E_i = 0, v_p = 0, E_p = 0) %>%
  #Calculate the probability of a symptomatic infection
  mutate(PR.sympt.inf = v * p_I * (1 - E_i) * s * (1 - E_s) + v_p * p_I * (1 - E_i) * s * (1 - E_p) + (1 - v - v_p) * p_I * s) %>%
  #Calculating baseline number of expected cases and declines under vaccination
  mutate(baseline_proportion = p_I*s,
         perc_decline = (baseline_proportion - PR.sympt.inf)/baseline_proportion*100,
         cases_averted_p1000 = baseline_proportion*1000 - PR.sympt.inf*1000,
         baseline_case_count_p1000 = baseline_proportion*1000,
         vaccinated_case_count_p1000 = PR.sympt.inf*1000) %>%
  #showing vaccinated proportion as a percentage for plotting
  mutate(v_perc = v*100)
  

p.cases <- dfc %>% 
  ggplot(aes(x = v, y = E_s, fill = cases_averted_p1000)) +
  geom_tile() +
  scale_fill_viridis_c(option = "turbo",
                       name = "Cases averted\n(per 1000\nindividuals\ntargeted)") +
  facet_grid(rows = vars(p_I), cols = vars(s)) +
  xlab(expression('Proportion fully vaccinated at time of disruptive event, v')) +
  ylab(expression('Vaccine efficacy against symptomatic infection, E'[S])) +
  theme_bw() +
  theme(axis.title = element_text(size = 16, family = "Arial"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        panel.grid = element_blank())
p.cases


#Plotting a subset for the main text
p.cases.trim <- dfc %>% 
  filter(s %in% c(0.4, 0.5)) %>%
  filter(p_I %in% c(0.3, 0.5, 0.7)) %>%
  ggplot(aes(x = v, y = E_s)) +
  geom_tile(aes(fill = cases_averted_p1000)) +
  geom_contour(aes(z = cases_averted_p1000), color = "black", breaks = c(100, 200, 300)) +
  scale_fill_viridis_c(option = "turbo",
                       name = "Cases averted\n(per 1000\nindividuals\ntargeted)") +
  facet_grid(rows = vars(p_I), cols = vars(s)) +
  xlab(expression('Proportion fully vaccinated at time of disruptive event, v')) +
  ylab(expression('Vaccine efficacy against symptomatic infection, E'[S])) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, family = "Arial"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 10), legend.title = element_text(size = 12),
        panel.grid = element_blank())
p.cases.trim

##################################################################
## Plotting cases averted for supplement
##################################################################

dfc.i.s <- tibble(expand.grid(v = (0:100)/100,
                            E = (0:100)/100,
                            I = c(0.1, 0.2, 0.3, 0.4),
                            #no cases expected if symptomatic infection = 0
                            s = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)))

dfc.s <- dfc.i.s %>%
  #Calculate the probability of a case (i.e., symptomatic infection) for each scenario
  mutate(PR.case = v*I*s*(1-E) + (1-v)*I*s) %>%
  #Calculating baseline number of expected cases and declines under vaccination
  mutate(baseline_proportion = I*s,
         perc_decline = (baseline_proportion - PR.case)/baseline_proportion*100,
         cases_averted_p1000 = baseline_proportion*1000 - PR.case*1000,
         baseline_case_count_p1000 = baseline_proportion*1000,
         vaccinated_case_count_p1000 = PR.case*1000)

p.cases.s <- dfc.s %>% 
  ggplot(aes(x = v, y = E, fill = cases_averted_p1000)) +
  geom_tile() +
  scale_fill_viridis_c(option = "turbo",
                       name = "Cases averted\n(per 1000 popn)") +
  facet_grid(rows = vars(I), cols = vars(s)) +
  labs(subtitle = "Probability infection is symptomatic (s) vs Probability of infection") +
  xlab(expression('Proportion fully vaccinated at time of disruptive event, v')) +
  ylab(expression('Vaccine efficacy against symptomatic infection, E'[S])) +
  theme_bw() +
  theme(axis.title = element_text(size = 16, family = "Arial"),
        strip.background = element_rect(fill = "white"),
        panel.grid = element_blank())
p.cases.s




p.cases.line.s <- dfc.s %>% filter(E == 0.75) %>% filter(v == 0.70) %>%
  ggplot(aes(x = s, y = cases_averted_p1000, color = I)) +
  geom_point() +
  scale_fill_viridis_c(option = "turbo",
                       name = "Cases averted\n(per 1000 popn)") +
  facet_grid(rows = vars(I), cols = vars(s)) +
  labs(subtitle = "Probability infection is symptomatic (s) vs Probability of infection") +
  xlab(expression('Proportion fully vaccinated at time of disruptive event, v')) +
  ylab(expression('Vaccine efficacy against symptomatic infection, E'[S])) +
  theme_bw() +
  theme(axis.title = element_text(size = 16, family = "Arial"),
        strip.background = element_rect(fill = "white"),
        panel.grid = element_blank())
p.cases.line.s


