library(tidyverse)
library(patchwork)



#Making an initial data frame for vaccination simulation
#using expand.grid() to write out all possible combinations
dfv.75.i <- tibble(expand.grid(v = 0.5, #fixed at 0.5 for this example
                            v.partial = (0:50)/100, #0-50% of unvaccinated individuals have some partial efficacy
                            s = 0.4, #fixed at 0.4 for this example
                            E = 0.75, #fixed at 0.75 for this example
                            E.partial = (0:75)/100)) #efficacy for those with partical vaccination completeness ranges from 0 to the same as fully vaccinated

dfv.75 <- dfv.75.i %>%
  #fixing the probability of infection at 0.365 for here
  mutate(I = 0.365) %>%
  #Calculate the probability of a symptomatic infection: fully vaccinated + partially vaccinated + unvaccinated
  mutate(PR.case = v*I*s*(1-E) + v.partial*I*s*(1-E.partial) + (1-v-v.partial)*I*s) %>%
  #Calculating baseline number of expected cases and declines under vaccination
  mutate(baseline_proportion = I*s,
         perc_decline = (baseline_proportion - PR.case)/baseline_proportion*100,
         cases_averted_p1000 = baseline_proportion*1000 - PR.case*1000,
         baseline_case_count_p1000 = baseline_proportion*1000,
         vaccinated_case_count_p1000 = PR.case*1000) %>%
  #showing vaccinated proportion as a percentage for plotting
  mutate(v_perc = v*100,
         v.partial_perc = v.partial*100)


#50% Coverage with R21
p.75 <- dfv.75 %>% 
  ggplot(aes(x = v.partial_perc, y = perc_decline, group = E.partial, color = E.partial)) +
  #Shading portion of the plot to higlight v and E combos with reduction > 50%
  annotate("rect", xmin=0, xmax=50, ymin=50, ymax=100, fill = "black", alpha = 0.07) +
  #Adding the geom_path
  geom_path(linewidth = 0.8, alpha = 0.8) +
  #50% reduction threshold
  geom_segment(aes(x=0, xend=50, y=50, yend=50), color="grey50", linetype="dashed", linewidth=0.3) +
  #Scales and themes
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  scale_y_reverse(limits = c(100, 0), breaks = seq(0, 100, by = 10)) +
  scale_color_viridis_c(option = "cividis",
                        name = "Avg. Efficacy for\nPartial Vaccinees",
                        end = 1,
                        limits = c(0, 0.75)) +
  labs(title = "R21 scenario (Vaccine efficacy = 0.75)") +
  xlab("Percent partially vaccinated at time of disruptive event (v)\n(in addition to 50% fully vaccinated)") +
  ylab(expression(paste("Percent reduction in symptomatic infections expected (", italic("vE"), ")"))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16, family = "Arial"),
        axis.text = element_text(size = 13, family = "Arial"))
p.75




#50% Coverage with a RTS,S

dfv.55.i <- tibble(expand.grid(v = 0.5,
                            v.partial = (0:50)/100,
                            s = 0.4,
                            E = 0.55, #fixing at 55% for RTS,S
                            E.partial = (0:55)/100))

dfv.55 <- dfv.55.i %>%
  #fixing the probability of infection at 0.365 for here
  mutate(I = 0.365) %>%
  #Calculate the probability of a symptomatic infection: fully vaccinated + partially vaccinated + unvaccinated
  mutate(PR.case = v*I*s*(1-E) + v.partial*I*s*(1-E.partial) + (1-v-v.partial)*I*s) %>%
  #Calculating baseline number of expected cases and declines under vaccination
  mutate(baseline_proportion = I*s,
         perc_decline = (baseline_proportion - PR.case)/baseline_proportion*100,
         cases_averted_p1000 = baseline_proportion*1000 - PR.case*1000,
         baseline_case_count_p1000 = baseline_proportion*1000,
         vaccinated_case_count_p1000 = PR.case*1000) %>%
  #showing vaccinated proportion as a percentage for plotting
  mutate(v_perc = v*100,
         v.partial_perc = v.partial*100)

p.55 <- dfv.55 %>% 
  ggplot(aes(x = v.partial_perc, y = perc_decline, group = E.partial, color = E.partial)) +
  #Shading portion of the plot to higlight v and E combos with reduction > 50%
  annotate("rect", xmin=0, xmax=50, ymin=50, ymax=100, fill = "black", alpha = 0.07) +
  #Adding the geom_path
  geom_path(linewidth = 0.8, alpha = 0.8) +
  #50% reduction threshold
  geom_segment(aes(x=0, xend=50, y=50, yend=50), color="grey50", linetype="dashed", linewidth=0.3) +
  #Scales and themes
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  scale_y_reverse(limits = c(100, 0), breaks = seq(0, 100, by = 10)) +
  scale_color_viridis_c(option = "cividis",
                        name = "Avg. Efficacy for\nPartial Vaccinees",
                        end = 1,
                        limits = c(0, 0.75)) +
  labs(title = "RTS,S scenario (Vaccine efficacy = 0.55)") +
  xlab("Percent partially vaccinated at time of disruptive event (v)\n(in addition to 50% fully vaccinated)") +
  ylab(expression(paste("Percent reduction in symptomatic infections expected (", italic("vE"), ")"))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16, family = "Arial"),
        axis.text = element_text(size = 13, family = "Arial"))
p.55


p.75 + p.55
