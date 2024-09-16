library(tidyverse)
library(patchwork)
library(ggrepel)


#Reading in data (corresponding to function call in determining FOI script)
dfi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_monthly.csv") %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10"))) %>%
  mutate(month = factor(month, levels = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"))) %>%
  mutate(month.num = month(ymd(paste0(2024, "-", month, "-", 1))))

p.foi <- dfi %>% 
  filter(day == 30) %>% 
  ggplot(aes(x = code.new, y = prob.inf, shape = age_cat, color = code.new, group = code.new)) +
  geom_path(color = "grey40", alpha = 0.4) +
  geom_point() +
  facet_wrap(vars(month), nrow = 2) +
  scale_color_viridis_d(option = "turbo",
                        name = "Site") +
  scale_shape_discrete(name = "Age Category") +
  xlab("Site") + ylab("Probability of infection") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        axis.text.x = element_text(angle = 90, family = "Andale Mono", hjust = 0.5, vjust = 0.5),
        axis.title = element_text(family = "Arial", size = 16),
        legend.position = "right")
p.foi


p.foi.month <- dfi %>% 
  filter(day == 30) %>% 
  ggplot(aes(x = month.num, y = prob.inf, color = age_cat)) +
  geom_path(aes(group = month.num), color = "grey40", alpha = 0.4) +
  geom_point() +
  facet_wrap(vars(code.new), nrow = 1) +
  scale_color_viridis_d(option = "turbo",
                        name = "Age Category") +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11)) +
  xlab("Month") + ylab("Probability of malaria infection") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = NULL),
        axis.text.x = element_text(angle = 90, family = "Andale Mono", 
                                   #size = 10,
                                   hjust = 1, vjust = 0.5),
        axis.title = element_text(family = "Arial", size = 16),
        legend.position = "bottom")
p.foi.month



#Looking up values
dfx <- dfi %>% filter(month == "MAY", day == "30", age == 13)
dfx <- dfi %>% filter(month == "OCT", day == "30", age == 13)

dfx <- dfi %>% filter(month == "MAY", site_code == "S1", age == 13)
