library(tidyverse)
library(patchwork)
library(ggrepel)


#Reading in data (corresponding to function call in determining FOI script)
dfi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_post_storm.csv") %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10"))) %>%
  mutate(age_cat = factor(age_cat, levels = c("Young children", "School aged children", "Teenagers", "Adults"))) %>%
  #Focusing on 60 days following the storms
  filter(day < 61)
  



p.line <- dfi %>% 
  filter(age %in% c(5, 13, 18, 27)) %>% 
  ggplot(aes(x = day, y = prob.inf, color = code.new, group = code.new)) +
  geom_path() +
  facet_grid(cols = vars(storm), rows = vars(age_cat)) +
  scale_color_viridis_d(option = "turbo",
                        name = "Site") +
  scale_x_continuous(breaks = c(0, 20, 40, 60)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6)) +
  xlab("Days after storm landfall") + ylab("Probability of infection") +
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NULL),
        axis.text = element_text(size = 10, family = "Andale Mono"),
        axis.text.x = element_text(family = "Andale Mono", hjust = 0.5, vjust = 0.5),
        axis.title = element_text(family = "Arial", size = 14),
        legend.position = "right")
p.line

p.heat.dark <- dfi %>% 
  ggplot(aes(x = day, y = age, fill = prob.inf)) +
  geom_tile(color = "black", linewidth = 0.02) +
  facet_grid(rows = vars(storm), cols = vars(code.new)) +
  scale_fill_viridis_c(option = "rocket",
                       name = "Probability of\ninfection") +
  scale_x_continuous(breaks = c(0, 20, 40, 60)) +
  xlab("Days after storm landfall") + ylab("Age") +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white", color = NULL),
        axis.text = element_text(size = 10, family = "Andale Mono"),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.title = element_text(family = "Arial", size = 14),
        legend.position = "right")
p.heat.dark

p.heat.dark + p.line + plot_layout(widths = c(10,3))


p.heat.light <- dfi %>% 
  ggplot(aes(x = day, y = age, fill = prob.inf)) +
  geom_tile() +
  facet_grid(rows = vars(storm), cols = vars(code.new)) +
  scale_fill_viridis_c(option = "inferno",
                       name = "Probability of\ninfection",
                       begin = 0.2) +
  scale_x_continuous(breaks = c(0, 20, 40, 60)) +
  xlab("Days after storm landfall") + ylab("Age") +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white", color = NULL),
        axis.text = element_text(size = 10, family = "Andale Mono"),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.title = element_text(family = "Arial", size = 14),
        legend.position = "right")
p.heat.light


p.heat.FRED <- dfi %>% 
  filter(storm == "CYCLONE FREDDY") %>%
  ggplot(aes(x = day, y = age, fill = prob.inf)) +
  geom_tile(color = "black", linewidth = 0.05) +
  facet_wrap(vars(code.new), nrow = 1) +
  scale_fill_viridis_c(option = "rocket",
                       name = "Probability of\ninfection") +
  scale_x_continuous(breaks = c(0, 20, 40, 60)) +
  xlab("Days after Cyclone Freddy landfall") + ylab("Age") +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white", color = NULL),
        axis.text = element_text(size = 14, family = "Arial"),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.title = element_text(family = "Arial", size = 20),
        strip.text = element_text(size = 16),
        legend.position = "right")
p.heat.FRED


#Looking up bounds
dfx <- dfi %>% filter(day == 60, age == 5)
dfx <- dfi %>% filter(day == 60, age == 13)

