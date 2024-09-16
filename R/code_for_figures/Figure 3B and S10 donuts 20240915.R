library(tidyverse)
library(geodata)
library(ggnewscale)
library(patchwork)
library(paletteer)

#########################################################################################################################
# Preparing data
#########################################################################################################################


#########################################################################################################################
# Data 1:Extracting temperature and precipitation data from worldclim
# using the geodata package to get monthly temperature and precipitation averages
# coordinates: approximate midpoint of Mananjary district coast

#writing a function to facilitate data extraction
f.extract.worldclim <- function(my.var, var.name.clean, my.lon, my.lat, xmin, xmax, ymin, ymax){
  #using worldclim_tile() to extract from the tile containing Mananjary
  raster.tile <- worldclim_tile(var = my.var, lon = my.lon, lat = my.lat, path = tempdir())
  my.names <- paste0(var.name.clean, 1:12)
  names(raster.tile) <- my.names
  #converting to data frame for plotting
  df1 <- as.data.frame(raster.tile, xy = TRUE, na.rm = TRUE)
  #trimming tile to approx. mid point of district
  df2 <- df1 %>% filter(x > xmin, x < xmax) %>% filter(y > ymin, y < ymax)
  df3 <- df2 %>% pivot_longer(!(x:y), names_to = "month.variable", values_to = "value") %>%
    #averaging across the rectangular area extracted
    group_by(month.variable) %>% summarize(mean = mean(value)) %>%
    #tidying data
    mutate(month.variable = factor(month.variable, levels = my.names)) %>% arrange(month.variable) %>%
    mutate(month = 1:12) %>% mutate(variable = my.var)
  return(df3)  
}


df.t.avg <- f.extract.worldclim(my.var         = "tavg",
                                var.name.clean = "t.avg.",
                                my.lon         = 48.3,
                                my.lat         = -21.3,
                                xmin           = 48.26118,
                                xmax           = 48.33878,
                                ymin           = -21.16840,
                                ymax           = -20.79915)

df.t.max <- f.extract.worldclim(my.var         = "tmax",
                                var.name.clean = "t.max.",
                                my.lon         = 48.3,
                                my.lat         = -21.3,
                                xmin           = 48.26118,
                                xmax           = 48.33878,
                                ymin           = -21.16840,
                                ymax           = -20.79915)

df.prec  <- f.extract.worldclim(my.var         = "prec",
                                var.name.clean = "prec.",
                                my.lon         = 48.3,
                                my.lat         = -21.3,
                                xmin           = 48.26118,
                                xmax           = 48.33878,
                                ymin           = -21.16840,
                                ymax           = -20.79915)


#Binding the dfs together
df.worldclim <- bind_rows(df.t.avg, df.t.max, df.prec)

#Test plotting climate variables
df.worldclim %>% filter(variable == "prec") %>%
  ggplot(aes(x = month, y = 10, fill = mean)) +
  #Ring 1: precip
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "prec", option = "mako") +
  #Ring 2: tavg
  new_scale_fill() +
  geom_tile(data = df.worldclim %>% filter(variable == "tavg"),
            aes(x = month, y = 8.9, fill = mean)) +
  scale_fill_viridis_c(name = "tavg", option = "magma", limits = c(20, 32)) +
  #Ring 3: tmax
  new_scale_fill() +
  geom_tile(data = df.worldclim %>% filter(variable == "tmax"),
            aes(x = month, y = 7.8, fill = mean)) +
  scale_fill_viridis_c(name = "tmax", option = "magma", limits = c(20, 32)) +
  
  ylim(0, 11) +
  coord_polar() +
  theme_void()


#########################################################################################################################
# Data 2:Cyclone landfall approach dates
# using file from output of a function call in the cyclone mapping (Figure 1) script

df.cycs <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_ibtracs_cyclones/processed/endemic_country_cyclone_hit_points.csv") %>%
  filter(COUNTRY == "Madagascar") %>%
  mutate(cyc.month = month(date), cyc.week = week(date), cyc.day = day(date))


#counts of storms by date
df.cycs.date_counts <- df.cycs %>%
  group_by(cyc.month, cyc.day) %>% summarize(cyc.n = n()) %>%
  mutate(date = ymd(paste0(2024, "-", cyc.month, "-", cyc.day)))

#Creating a data frame with all dates in a calendar year as rows (using 2024 as an example for convenience due to leap years)
dfm.day <- tibble(date = seq(ymd('2024-01-01'),ymd('2024-12-31'), by = 'days'),
              month = month(date),
              day = day(date)) %>%
  left_join(df.cycs.date_counts, by = join_by(date)) %>%
  #Setting dates with zero storms to 0
  mutate(cyc.n = ifelse(is.na(cyc.n), 0, cyc.n)) %>%
  mutate(cyc.n.fac = factor(cyc.n, levels = c("0", "1", "2", "3", "4")))

dfm.day %>% 
  ggplot(aes(x = date, y = 10, fill = cyc.n.fac)) +
  geom_tile(color = "white", linewidth = 0.05) +
  geom_point(aes(x = date, y = 10.6), color = "black", size = 0.1) +
  scale_fill_viridis_d(name = "Cyclone Count", 
                       option = "mako") +
  ylim(0, 11.5) +
  coord_polar() +
  theme_void()

#########################################################################################################################
# Joining cyclone and climate data

#tidying climate data
df.world.clim.c <- df.worldclim %>% 
  #keeping avg temp only
  filter(variable != "tmax") %>%
  #pivoting wider
  dplyr::select(-month.variable) %>% pivot_wider(names_from = "variable", values_from = "mean")

dfm.i <- dfm.day %>% 
  mutate(tavg = case_when(
    month ==  1 ~ df.world.clim.c$tavg[1],
    month ==  2 ~ df.world.clim.c$tavg[2],
    month ==  3 ~ df.world.clim.c$tavg[3],
    month ==  4 ~ df.world.clim.c$tavg[4],
    month ==  5 ~ df.world.clim.c$tavg[5],
    month ==  6 ~ df.world.clim.c$tavg[6],
    month ==  7 ~ df.world.clim.c$tavg[7],
    month ==  8 ~ df.world.clim.c$tavg[8],
    month ==  9 ~ df.world.clim.c$tavg[9],
    month == 10 ~ df.world.clim.c$tavg[10],
    month == 11 ~ df.world.clim.c$tavg[11],
    month == 12 ~ df.world.clim.c$tavg[12])) %>%
  mutate(prec = case_when(
    month ==  1 ~ df.world.clim.c$prec[1],
    month ==  2 ~ df.world.clim.c$prec[2],
    month ==  3 ~ df.world.clim.c$prec[3],
    month ==  4 ~ df.world.clim.c$prec[4],
    month ==  5 ~ df.world.clim.c$prec[5],
    month ==  6 ~ df.world.clim.c$prec[6],
    month ==  7 ~ df.world.clim.c$prec[7],
    month ==  8 ~ df.world.clim.c$prec[8],
    month ==  9 ~ df.world.clim.c$prec[9],
    month == 10 ~ df.world.clim.c$prec[10],
    month == 11 ~ df.world.clim.c$prec[11],
    month == 12 ~ df.world.clim.c$prec[12]))

#Test plotting
dfm.i %>% 
  ggplot(aes(x = date, y = 10, fill = cyc.n.fac)) +
  #Ring 1: cyclone day counts
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(name = "Historical Cyclone\nCount", values = c("white", "#F5B6C8", "#DE89B0", "#C868A1", "#9C4286")) +
  #Ring 2: precip
  new_scale_fill() +
  geom_tile(aes(x = date, y = 8.9, fill = prec), linewidth = 0.1) +
  scale_fill_viridis_c(name = "Precip (mm)", option = "mako", begin = 0.1) +
  #Ring 3: tavg
  new_scale_fill() +
  geom_tile(aes(x = date, y = 7.8, fill = tavg), linewidth = 0.1) +
  scale_fill_viridis_c(name = "Avg Temp (C)", option = "magma", begin = 0.1) +
  #Outer labels:
  #new_scale_fill() +
  #geom_point(aes(x = date, y = 10.7), color = "black", size = 0.1) +
  #Theme and coord params
  ylim(0, 11) +
  coord_polar() +
  theme_void()

#########################################################################################################################
# Adding malaria data
df.foi <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_foi/foi_bookends 20231023.csv") %>%
  dplyr::select(date, S6.14) %>% rename(prob.inf = S6.14) %>%
  #note 2024 was a leap year so will need to account for Feb 29
  mutate(date.2024 = ymd(paste0(2024, "-", month(date), "-", day(date)))) %>% dplyr::select(-date)

dfm <- dfm.i %>% left_join(df.foi, by = join_by(date == date.2024)) %>%
  #interpolating February 29th simply for visualization: using midpoint between day before and day after
  mutate(prob.inf = ifelse(date == "2024-02-29", (0.2612645 + 0.2633793)/2, prob.inf)) %>%
  #adding a month label
  mutate(month.label = month(date, label = TRUE, abbr = TRUE))

#Making a vector of the date halfway through a month to facilitate labeling
v.mid.months <- ymd(c("2024-01-15",
                      "2024-02-14",
                      "2024-03-15",
                      "2024-04-15",
                      "2024-05-15",
                      "2024-06-15",
                      "2024-07-15",
                      "2024-08-15",
                      "2024-09-15",
                      "2024-10-15",
                      "2024-11-15",
                      "2024-12-15"))

#Making a vector of the date halfway through a month to facilitate labeling
v.months.bounds <- ymd(c("2024-01-01",
                         "2024-02-01",
                         "2024-03-01",
                         "2024-04-01",
                         "2024-05-01",
                         "2024-06-01",
                         "2024-07-01",
                         "2024-08-01",
                         "2024-09-01",
                         "2024-10-01",
                         "2024-11-01",
                         "2024-12-01"))


#########################################################################################################################
# Figure 3B: Cyclones-Precip-Temp-Malaria (for an exemplar site)
#########################################################################################################################
dfm %>% 
  ggplot(aes(x = date, y = 10, fill = cyc.n.fac)) +
  #Ring 1: cyclone day counts
  geom_tile(color = "white", height = 1.2, width = 0.8) +
  scale_fill_manual(name = "Historical Cyclone\nCount", values = c("white", "#F5B6C8", "#DE89B0", "#C868A1", "#9C4286"),
                    guide = guide_legend(order = 1)) +
  #Ring 2: precip
  new_scale_fill() +
  geom_tile(aes(x = date, y = 9.1, fill = prec), height = 0.5, width = 1.1) +
  scale_fill_paletteer_c("grDevices::Teal", name = "Precip (mm)", direction = -1, guide = guide_colourbar(order = 2)) +
  #Ring 3: tavg
  new_scale_fill() +
  geom_tile(aes(x = date, y = 8.5, fill = tavg), height = 0.5, width = 1.1) +
  scale_fill_paletteer_c("grDevices::Heat", name = "Avg Temp (C)", direction = -1, guide = guide_colourbar(order = 3)) +
  #Ring 4: malaria
  new_scale_fill() +
  geom_tile(aes(x = date, y = 7.35, fill = prob.inf), width = 1.2, height = 1.6) +
  #Color scale
  scale_fill_paletteer_c("grDevices::Blue-Red 3", name = "Cumulative probability\nof malaria infection\n(per 31 days)") +
  

  #Outer labels:
  new_scale_fill() +
  #Dots for Cyclones Freddy and Batsirai
  geom_point(data = dfm %>% filter(date %in% c(ymd("2024-02-05"), ymd("2024-2-21"))),
             aes(x = date, y = 10.85),
             color = "#C868A1") +
  #Ticks for each day along 'x' axis
  geom_segment(aes(x = date, xend = date, y = 6.35, yend = 6.5),
               color = "black", linewidth = 0.25, alpha = 1, linetype = "solid") +
  #Month boundaries
  geom_segment(data = dfm %>% filter(date %in% v.months.bounds),
               aes(x = date, xend = date, y = 5, yend = 6),
               color = "black", linewidth = 0.3, alpha = 0.5, linetype = "solid") +
  #Month labels
  geom_text(data = dfm %>% filter(date %in% v.mid.months),
            aes(x = date, y = 5.48, label = month.label),
            color = "black", alpha = 0.8, size = 3.5) +
  #Theme and coord params
  ylim(0, 11) +
  coord_polar() +
  theme_void() +
  theme(legend.key = element_rect(color = "black"),
        #legend.key.size = unit(1.1, 'cm'),
        #legend.title = element_text(size = 14)
        )


#########################################################################################################################
# Figure S10: Seasonality for all sites
#########################################################################################################################

df.foi.all <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_foi/foi_bookends 20231023.csv") %>%
  pivot_longer(!date, names_to = "foi.site", values_to = "prob.inf") %>%
  mutate(site_code = substr(foi.site, 1, 2)) %>% 
  mutate(age = as.integer(substr(foi.site, 4, 5))) %>%
  #updating site codes
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
  #Ordering by average infection rate
  mutate(code.new = factor(code.new, levels = c("MNJ.06",
                                                "MNJ.04",
                                                "MNJ.03",
                                                "MNJ.05",
                                                "MNJ.02",
                                                "MNJ.07",
                                                "MNJ.08",
                                                "MNJ.09",
                                                "MNJ.01",
                                                "MNJ.10"))) %>%
  mutate(month = month(date)) %>%
  mutate(month.label = month(date, label = TRUE, abbr = TRUE))

#Making a vector of the date halfway through a month to facilitate labeling
v.mid.months.2023 <- ymd(c("2023-01-15",
                           "2023-02-14",
                           "2023-03-15",
                           "2023-04-15",
                           "2023-05-15",
                           "2023-06-15",
                           "2023-07-15",
                           "2023-08-15",
                           "2023-09-15",
                           "2023-10-15",
                           "2023-11-15",
                           "2023-12-15"))

#Making a vector of the date halfway through a month to facilitate labeling
v.months.bounds.2023 <- ymd(c("2023-01-01",
                              "2023-02-01",
                              "2023-03-01",
                              "2023-04-01",
                              "2023-05-01",
                              "2023-06-01",
                              "2023-07-01",
                              "2023-08-01",
                              "2023-09-01",
                              "2023-10-01",
                              "2023-11-01",
                              "2023-12-01"))
#Donut 1: Young children
p.d.05 <- df.foi.all %>% filter(age == 5) %>%
  ggplot(aes(x = date, y = 20, fill = prob.inf)) +
  #Site 1
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.06") %>% filter(age == 5),  width = 1.2, height = 1) +
  #Sites 2-10
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.04") %>% filter(age == 5), aes(x = date, y = 19), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.03") %>% filter(age == 5), aes(x = date, y = 18), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.05") %>% filter(age == 5), aes(x = date, y = 17), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.02") %>% filter(age == 5), aes(x = date, y = 16), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.07") %>% filter(age == 5), aes(x = date, y = 15), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.08") %>% filter(age == 5), aes(x = date, y = 14), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.09") %>% filter(age == 5), aes(x = date, y = 13), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.01") %>% filter(age == 5), aes(x = date, y = 12), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.10") %>% filter(age == 5), aes(x = date, y = 11), width = 1.2, height = 1) +
  
  #Color scale
  scale_fill_paletteer_c("grDevices::Blue-Red 3", name = "Cumulative probability\nof malaria infection\n(per 31 days)", limits = c(0, max(df.foi.all$prob.inf))) +
  
  #Inner labels:
  new_scale_fill() +
  geom_text(data = df.foi.all %>% filter(age == 5) %>% filter(date == ymd("2023-01-01")),
            aes(x = date, 
                y = 20:11, 
                label = c("MNJ.06","MNJ.04","MNJ.03","MNJ.05","MNJ.02","MNJ.07","MNJ.08","MNJ.09","MNJ.01","MNJ.10")),
            size = 2.5, color = "white") +
  
  #Outer labels:
  new_scale_fill() +
  #Ticks for each day along 'x' axis
  geom_segment(data = df.foi.all %>% filter(code.new == unique(dfm.all$code.new)[1])  %>% filter(age == 5),
               aes(x = date, xend = date, y = 20.6, yend = 20.9),
               color = "black", linewidth = 0.25, alpha = 1, linetype = "solid") +
  #Month boundaries
  geom_segment(data = df.foi.all %>% filter(date %in% v.months.bounds.2023) %>% filter(code.new == unique(dfm.all$code.new)[1]) %>% filter(age == 5),
               aes(x = date, xend = date, y = 21.2, yend = 22),
               color = "black", linewidth = 0.3, alpha = 0.5, linetype = "solid") +
  #Month labels
  geom_text(data = df.foi.all %>% filter(date %in% v.mid.months.2023) %>% filter(code.new == unique(dfm.all$code.new)[1]) %>% filter(age == 5),
            aes(x = date, y = 21.8, label = month),
            color = "black") +
  #Theme and coord params
  ylim(0, 22) +
  coord_polar() +
  theme_void() +
  theme(legend.position = "bottom")
p.d.05

#Donut 2: School aged children
p.d.14 <- df.foi.all %>% filter(age == 14) %>%
  ggplot(aes(x = date, y = 20, fill = prob.inf)) +
  #Site 1
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.06") %>% filter(age == 14),  width = 1.2, height = 1) +
  #Sites 2-10
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.04") %>% filter(age == 14), aes(x = date, y = 19), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.03") %>% filter(age == 14), aes(x = date, y = 18), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.05") %>% filter(age == 14), aes(x = date, y = 17), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.02") %>% filter(age == 14), aes(x = date, y = 16), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.07") %>% filter(age == 14), aes(x = date, y = 15), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.08") %>% filter(age == 14), aes(x = date, y = 14), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.09") %>% filter(age == 14), aes(x = date, y = 13), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.01") %>% filter(age == 14), aes(x = date, y = 12), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.10") %>% filter(age == 14), aes(x = date, y = 11), width = 1.2, height = 1) +
  
  #Color scale
  scale_fill_paletteer_c("grDevices::Blue-Red 3", name = "Cumulative probability\nof malaria infection\n(per 31 days)", limits = c(0, max(df.foi.all$prob.inf))) +
  
  #Outer labels:
  new_scale_fill() +
  #Ticks for each day along 'x' axis
  geom_segment(data = df.foi.all %>% filter(code.new == unique(dfm.all$code.new)[1])  %>% filter(age == 14),
               aes(x = date, xend = date, y = 20.6, yend = 20.9),
               color = "black", linewidth = 0.25, alpha = 1, linetype = "solid") +
  #Month boundaries
  geom_segment(data = df.foi.all %>% filter(date %in% v.months.bounds.2023) %>% filter(code.new == unique(dfm.all$code.new)[1]) %>% filter(age == 14),
               aes(x = date, xend = date, y = 21.2, yend = 22),
               color = "black", linewidth = 0.3, alpha = 0.5, linetype = "solid") +
  #Month labels
  geom_text(data = df.foi.all %>% filter(date %in% v.mid.months.2023) %>% filter(code.new == unique(dfm.all$code.new)[1]) %>% filter(age == 14),
            aes(x = date, y = 21.8, label = month),
            color = "black") +
  #Theme and coord params
  ylim(0, 22) +
  coord_polar() +
  theme_void() +
  theme(legend.position = "bottom")
p.d.14
#Donut 3: Adults
p.d.30 <- df.foi.all %>% filter(age == 30) %>%
  ggplot(aes(x = date, y = 20, fill = prob.inf)) +
  #Site 1
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.06") %>% filter(age == 30),  width = 1.2, height = 1) +
  #Sites 2-10
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.04") %>% filter(age == 30), aes(x = date, y = 19), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.03") %>% filter(age == 30), aes(x = date, y = 18), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.05") %>% filter(age == 30), aes(x = date, y = 17), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.02") %>% filter(age == 30), aes(x = date, y = 16), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.07") %>% filter(age == 30), aes(x = date, y = 15), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.08") %>% filter(age == 30), aes(x = date, y = 14), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.09") %>% filter(age == 30), aes(x = date, y = 13), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.01") %>% filter(age == 30), aes(x = date, y = 12), width = 1.2, height = 1) +
  geom_tile(data = df.foi.all %>% filter(code.new == "MNJ.10") %>% filter(age == 30), aes(x = date, y = 11), width = 1.2, height = 1) +
  
  #Color scale
  scale_fill_paletteer_c("grDevices::Blue-Red 3", name = "Cumulative probability\nof malaria infection\n(per 31 days)", limits = c(0, max(df.foi.all$prob.inf))) +
  
  #Outer labels:
  new_scale_fill() +
  #Ticks for each day along 'x' axis
  geom_segment(data = df.foi.all %>% filter(code.new == unique(dfm.all$code.new)[1])  %>% filter(age == 30),
               aes(x = date, xend = date, y = 20.6, yend = 20.9),
               color = "black", linewidth = 0.25, alpha = 1, linetype = "solid") +
  #Month boundaries
  geom_segment(data = df.foi.all %>% filter(date %in% v.months.bounds.2023) %>% filter(code.new == unique(dfm.all$code.new)[1]) %>% filter(age == 30),
               aes(x = date, xend = date, y = 21.2, yend = 22),
               color = "black", linewidth = 0.3, alpha = 0.5, linetype = "solid") +
  #Month labels
  geom_text(data = df.foi.all %>% filter(date %in% v.mid.months.2023) %>% filter(code.new == unique(dfm.all$code.new)[1]) %>% filter(age == 30),
            aes(x = date, y = 21.8, label = month),
            color = "black") +
  #Theme and coord params
  ylim(0, 22) +
  coord_polar() +
  theme_void() +
  theme(legend.position = "bottom")
p.d.30

p.d.05 + p.d.14 + p.d.30











