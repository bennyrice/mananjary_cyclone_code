library(tidyverse)
library(patchwork)
library(ggrepel)


#Reading in data (corresponding to function call in determining FOI script)
dfi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/R/modeling/output/foi_by_age.csv") %>%
  mutate(code.new = factor(code.new, levels = c("MNJ.06", "MNJ.04", "MNJ.03", "MNJ.05", "MNJ.02", "MNJ.07", "MNJ.08", "MNJ.09", "MNJ.01", "MNJ.10")))


p.foi.age <- dfi %>% 
  #For visualization, showing the month of January for all sites
  #Age range trimmed to less than 75 yos due to paucity of enrollees over age 75
  filter(month == "JAN") %>% filter(age < 75) %>%
  ggplot(aes(x = age, y = prob.inf, color = code.new, group = code.new)) +
  geom_point() +
  geom_path() +
  scale_color_viridis_d(option = "turbo", name = "Site Code") +
  xlab("Age (years)") + ylab("Force of infection (FOI)\n(per day)") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "Andale Mono", size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"))
p.foi.age