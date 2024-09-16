library(tidyverse)
library(patchwork)
library(ggrepel)
library(sf)
library(biscale)
library(raster)
library(scales)


# 1a  Global map of malaria incidence surface + cyclone tracks
# 1ai Year and intensity for the 1705 recorded storms impacting malaria endemic countries since 1980

# 1b  Bivariate color map: Cyclone exposure vs malaria incidence
# 1bi Scatter plot of malaria incidence vs cyclone activity at national level

# 1c  Cyclone tracks in Africa
# 1ci Cyclone impacts by African country

# 1d  Cyclones in SE Madagascar and Mananjary district

###########################################################################################
## Reading in data up front
###########################################################################################

###########################################################################################
#Data 1: Reading in a list of countries and ISO codes for later joining
df.iso <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_iso_codes/iso_codes.csv")

###########################################################################################
#Data 2: Reading in a global country shape file

#downloaded from https://hub.arcgis.com/datasets/esri::world-countries-generalized/explore
countries <- st_read("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_global_shape_files/global_shape_file/World_Countries_Generalized.shp")
#*note that Taiwan and Hong Kong are included in the geometry for China it appears
#Cleaning up replicate ISO2 codes
countries <- countries %>% 
  mutate(ISO = case_when(
    COUNTRY == "Namibia" ~ "NAm", #otherwise can appear as a no data 'NA' later
    
    COUNTRY == "Bonaire"         ~ "BQb",
    COUNTRY == "Saba"            ~ "BQs",
    COUNTRY == "Saint Eustatius" ~ "BQe",
    
    COUNTRY == "Canarias" ~ "ESc",
    
    COUNTRY == "French Southern Territories" ~ "TFs",
    COUNTRY == "Glorioso Islands"            ~ "TFg",
    COUNTRY == "Juan De Nova Island"         ~ "TFj",
    
    COUNTRY == "Azores"   ~ "PTa",
    COUNTRY == "Madeira"  ~ "PTm",
    COUNTRY == "Portugal" ~ "PT",
    
    .default = as.character(ISO)))

#joining with ISO codes and standard names
df.countries <- df.iso %>% full_join(countries, by = join_by(ISO2 == ISO)) %>%
  #Dropping islands etc already included in a country geometry
  filter(!is.na(SHAPE_Leng))

###########################################################################################
#Data 3: Reading in cyclone affected data from Jing et al 2023
#PUBMED Link: https://pubmed.ncbi.nlm.nih.gov/38122822/
#Country ISO codes added to Supplementary table 1 from the paper for later joining
df.jing <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_cyclone_exposure/Jing et al 2024 Supplement/supplementary_table1c.csv")

###########################################################################################
#Data 4: Reading in malaria incidence data from Malaria Atlas Project 2022 raster

#Reading in, converting to data frame then plot using geom_raster
#Read in using raster()
r.mal.ap <- raster::raster("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_MAP_incidence_raster/2024_GBD2023_Global_Pf_Incidence_Rate_2000/2024_GBD2023_Global_Pf_Incidence_Rate_2022.tif")
#Convert to data frame
df.mal.ap <- as.data.frame(r.mal.ap, xy = TRUE) %>% 
  #cleaning up column names
  rename(inc_2022 = X2024_GBD2023_Global_Pf_Incidence_Rate_2022_1) %>% 
  #scale to incidence per 1000 popn
  mutate(inc_2022 = 1000*inc_2022) %>%
  #Convert 0s to NAs
  mutate(inc_2022 = ifelse(inc_2022 == 0, NA, inc_2022)) %>%
  #Rescaling near elimination settings to 0.1 per MAP
  # mutate(inc_2020 = ifelse(inc_2020 <= 1, 1, inc_2020))
  #Rescaling low transmission settings to 1 for downstream plotting on a legible log scale
  mutate(inc_2022 = ifelse(inc_2022 <= 1, 1, inc_2022))

###########################################################################################
#Data 5: Reading in cyclone tracks data from IBTRACs
dfa.i <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_ibtracs_cyclones/ibtracs.since1980.list.v04r01.csv")

cols_to_keep.a <- c("SEASON", 
                    "SID",
                    "NAME", 
                    "ISO_TIME", 
                    "LAT", "LON", 
                    "DIST2LAND", #in km
                    "WMO_WIND", "USA_WIND" #in knots (note sometimes WMO or USA don't have wind data)
)

#Cleaning up the storm tracks data
dfa <- dfa.i %>% 
  #subset to columns of interest
  dplyr::select(all_of(cols_to_keep.a)) %>%
  #drop the dummy row that is comments on units from the ibtracs file
  filter(!is.na(NAME)) %>% 
  #convert to numeric where needed
  mutate_at(c('SEASON', 'LAT', 'LON', 'DIST2LAND', 'WMO_WIND', 'USA_WIND'), as.numeric) %>%
  #Creating a unique name for each storm
  mutate(name.yr = paste0(SEASON, "_", SID, "_", NAME)) %>%
  #Dealing with lon values outside the expected range -180 to 180
  mutate(LON = ifelse(LON > 180, (360 %% LON) * -1, LON)) %>%
  #Dropping storms without a track (ie only crossed threshold for inclusion for 1 data point)
  group_by(name.yr) %>% mutate(n_points = n()) %>% ungroup() %>% filter(n_points > 1) %>%
  #Because of the way wrapping around the international date line is handled, need to offset from 180 slightly
  mutate(LON = ifelse(LON == 180, 179.999, LON)) %>%
  #converting wind speed to kmh
  mutate(wind.speed = USA_WIND*1.852)

###########################################################################################
#Data 6: Reading in Madagascar District Shape Files

#Read in District (ADM2) boundaries from shape file
mg.adm2 <- st_read("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_madagascar_admin/mdg_adm_bngrc_ocha_20181031_shp/mdg_admbnda_adm2_BNGRC_OCHA_20181031.shp")
#Read in Country (ADM0) boundaries for Madagascar
mg.adm0 <- st_read("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_madagascar_admin/mdg_adm_bngrc_ocha_20181031_shp/mdg_admbnda_adm0_BNGRC_OCHA_20181031.shp")

###########################################################################################
#Data 6: Reading in national malaria incidence estimates from WHO
#Downloaded 20240802 from WHO Global Health Observatory (GHO)
#   incidence (inc):  "Estimated malaria incidence (per 1000 population at risk)"
#   cases (cases):    "Estimated number of malaria cases"

df.inc.who <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_WHO_malaria_incidence_data/WHO incidence data 20240802.csv") %>%
  #Trimming to most recent data
  filter(IsLatestYear == TRUE) %>%
  #Tidying relevant variables
  dplyr::select(Location, SpatialDimValueCode, Period, FactValueNumeric) %>%
  rename(country = Location, ISO3 = SpatialDimValueCode, inc.year = Period, inc = FactValueNumeric)

df.cases.who <- readr::read_csv("/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_mapping_WHO_malaria_incidence_data/WHO case count data 20240802.csv") %>%
  #Trimming to most recent data
  filter(IsLatestYear == TRUE) %>%
  #Tidying relevant variables
  dplyr::select(Location, SpatialDimValueCode, Period, FactValueNumeric) %>%
  rename(country = Location, ISO3 = SpatialDimValueCode, cases.year = Period, cases = FactValueNumeric)

#Joining WHO incidence data with global countries shape file for ease of plotting later
df.countries.inc <- df.countries %>% left_join(df.inc.who, by = join_by(ISO3)) %>% 
  #dropping redundant column
  dplyr::select(-country) %>%
  #cleaning up incidence column
  mutate(inc = ifelse(is.na(inc), 0, inc)) %>%
  #Adding case count data
  left_join(df.cases.who, by = join_by(ISO3)) %>%
  #cleaning up case count column
  mutate(cases = ifelse(is.na(cases), 0, cases))


###########################################################################################
## Figure 1a ## Global map of malaria incidence surface + cyclone tracks
###########################################################################################

#Calculating minimum distance to land for a storm
dfa.land.i <- dfa %>% group_by(SEASON, name.yr) %>% summarize(min.DIST2LAND = min(DIST2LAND)) 
#Exploratory plot of proximity to land for all storms
dfa.land.i %>%
  ggplot(aes(x = SEASON, y = min.DIST2LAND, color = name.yr)) +
  geom_jitter(alpha = 0.5) + 
  scale_color_viridis_d(option = "mako") +
  theme_bw() + theme(legend.position = "none")
#Exploratory plot of proximity to land for storms within 112 km of land
dfa.land.i %>% filter(min.DIST2LAND < 112) %>%
  ggplot(aes(x = SEASON, y = min.DIST2LAND, color = name.yr)) +
  geom_jitter(alpha = 0.5) + 
  scale_color_viridis_d(option = "mako") +
  theme_bw() + theme(legend.position = "none")

#Converting to LINESTRING / MULTISTRING to account for date line
#First, convert data frame with lat and lon to a sf object of coordinate points and set projection
sfa.i <- st_as_sf(dfa, coords = c("LON", "LAT")) %>% st_set_crs(4326)
#Second, group sets of points into a line
sfa <- sfa.i %>% group_by(name.yr) %>% summarise(do_union = FALSE) %>% st_cast("LINESTRING")
#Third, break lines into multilines to wrap around intl date line
sfa.t <- st_wrap_dateline(sfa, options = "WRAPDATELINE=YES", quiet = TRUE)


#### Plotting malaria incidence with storm tracks
p.m.c.inc_map <- ggplot() + 
  #Cyclone layer: Adding cyclone tracks
  #(*eventually subsetting or counting those that contact a malaria endemic locale)
  #(*3 categories of storms, those that don't make land, make land (non-endemic), make land (malaria endemic))
  #STORM TRACKS: LINE
  geom_sf(data = sfa.t, color = "#358AAA", alpha = 0.1, linewidth = 0.5) +
  #Base layer: gobal map from country shape files, dropping Antartica
  geom_sf(data = countries %>% filter(ISO != "AQ"),
          fill = "black", alpha = 0.92, color = "black", size = 0.01, show.legend = FALSE) +
  #Malaria layer: Adding layer for incidence data from MAP using geom_raster
  geom_raster(data = df.mal.ap, aes(x = x, y = y, fill = inc_2022)) +
  #Using Viridis:Rocket, upping begin to remove some black, and log scale to deal with many lows
  scale_fill_viridis_c(option = "rocket",
                       #na.value set to transparent to not show patches of raster with no malaria
                       na.value = "transparent",
                       begin = 0.05,
                       end = 0.95,
                       transform = "log2",
                       breaks = c(1, 10, 100, 500),
                       name = "Malaria\nincidence\n(per 1,000 popn)") +
  #Set coords
  coord_sf() +
  theme(legend.position = "right",
        #legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))
p.m.c.inc_map





###########################################################################
## Figure 1B ## Number of storms hitting malaria endemic countries
###########################################################################

#Identifying malaria endemic countries
v.endemic_countries <- df.countries.inc$ISO3[df.countries.inc$inc > 0]

df.endemic <- df.countries.inc %>% filter(inc > 0)

sf.endemic <- st_as_sf(x = df.endemic, sf_column_name = "geometry") %>% st_set_crs(4326)

#Making a 60nm = 111km buffer around malaria endemic countries
sf.END.buff <- st_buffer(sf.endemic, dist = 111*1000)

#Using st_intersection to find storms tracks that pass within buffer of malaria endemic countries
int.sf.END <- st_intersection(sfa.t, sf.END.buff)

v.endemic.hits <- unique(int.sf.END$name.yr)

df.endemic.hits <- dfa %>% filter(name.yr %in% v.endemic.hits) %>%
  #dropping the n column for clarity (will be counting data points on subsequent groups)
  dplyr::select(-n_points) %>%
  #filtering for storms within 60 nm of land
  filter(DIST2LAND < 112) %>% 
  #convert date to date format with lubridate
  mutate(date = ymd(substr(ISO_TIME, 1, 10))) %>%
  #Finding the max wind speed on storm approach (within 60nm of land)
  group_by(name.yr) %>% 
  mutate(n = n(), n_NAs = sum(is.na(wind.speed)), n_wind.points = n() - n_NAs) %>%
  #For storms with wind speed data, find point with max wind
  mutate(max.wind.speed = ifelse(n_wind.points > 0, max(wind.speed, na.rm = TRUE), NA)) %>% ungroup() %>%
  filter((wind.speed == max.wind.speed) %>% replace_na(TRUE)) %>%
  #When multiple points have the same max wind, choose the one closest to land and earliest date
  group_by(name.yr) %>%
  mutate(min.DIST2LAND = min(DIST2LAND)) %>% ungroup() %>%
  filter(DIST2LAND == min.DIST2LAND) %>%
  group_by(name.yr) %>%
  mutate(min.date = min(date)) %>% ungroup() %>%
  filter(date == min.date) %>%
  #When max wind was the same at multiple points on earliest date, arbitrarily keep one data point
  group_by(name.yr) %>%
  arrange(ISO_TIME) %>%
  filter(row_number() == 1) %>%
  mutate(n.points = n()) %>% ungroup() %>%
  #Adding the year of the hit
  mutate(year = year(date)) %>%
  #Adding storm severity category
  mutate(storm.cat = case_when(
    max.wind.speed < 62                              ~ "TD",
    max.wind.speed > 62    & max.wind.speed <= 118.5 ~ "TS",
    max.wind.speed > 118.5 & max.wind.speed <= 153   ~ "H1",
    max.wind.speed > 153   & max.wind.speed <= 177   ~ "H2",
    max.wind.speed > 177   & max.wind.speed <= 208   ~ "H3",
    max.wind.speed > 208   & max.wind.speed <= 251   ~ "H4",
    max.wind.speed > 251                             ~ "H5",
    is.na(max.wind.speed) ~ "ND")) %>%
  #setting factor levels
  mutate(storm.cat = factor(storm.cat, levels = c("ND", "TD", "TS", "H1", "H2", "H3", "H4", "H5")))


#Plot 1: Bar charts by storm category over time
p.endemic.hits.bar.by_cat <- df.endemic.hits %>% 
  group_by(year, storm.cat) %>%
  summarize(n = n()) %>% ungroup() %>%
  #drop 2024 as data is not yet complete
  #filter(year < 2024) %>%
  ggplot(aes(x = year, y = n, fill = storm.cat)) +
  geom_bar(position = "stack", stat = "identity") + 
  facet_wrap(vars(storm.cat), nrow = 8, strip.position = "right") +
  scale_fill_viridis_d(option = "magma",
                       name = "Storm\nCategory",
                       direction = 1,
                       end = 0.90) +
  xlab("Year") +
  ylab("Number of storms per year among currently malaria endemic countries") +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", angle = 0, hjust = 0.5, vjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(family = "Arial", size = 14),
        strip.background = element_rect(fill = "white"))
p.endemic.hits.bar.by_cat


#Plot 2: Box and whisker plot by storm category
p.endemic.hits.by_cat <- df.endemic.hits %>% group_by(storm.cat, year) %>% summarize(n = n()) %>%
  ggplot(aes(x = storm.cat, y = n, color = storm.cat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.18, height = 0.2, alpha = 0.9) +
  scale_color_viridis_d(option = "mako",
                        name = "Storm\nCategory",
                        end = 0.88) +
  xlab("Storm Category") + 
  ylab("Number of storms per year among currently malaria endemic countries\n(since 1980)") +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", angle = 0, hjust = 0.5, vjust = 0.5),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.title = element_text(family = "Arial", size = 18),
        axis.text = element_text(family = "Arial", size = 14))
p.endemic.hits.by_cat

#Plot 3a: Storms since 1980 colored by storm category
p.endemic.hits.bar.by_cat <- df.endemic.hits %>% 
  group_by(year, storm.cat) %>%
  summarize(n = n()) %>% ungroup() %>%
  #drop 2024 as data is not yet complete
  #filter(year < 2024) %>%
  ggplot(aes(x = year, y = n, fill = storm.cat)) +
  geom_bar(position = "stack", stat = "identity") + 
  #facet_wrap(vars(storm.cat), nrow = 8, strip.position = "right") +
  scale_fill_viridis_d(option = "mako",
                       name = "Storm\nCategory",
                       end = 0.88) +
  xlab("Year") +
  ylab("Number of storms per year among currently malaria endemic countries") +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", angle = 0, hjust = 0.5, vjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(family = "Arial", size = 18),
        axis.text = element_text(family = "Arial", size = 14),
        strip.background = element_rect(fill = "white"))
p.endemic.hits.bar.by_cat


p.endemic.hits.by_cat + p.endemic.hits.bar.by_cat


###########################################################################################
## Supplemental Figure S1A ## Making a bivariate color map of cyclones and malaria
###########################################################################################

#malaria incidence data with shape files
# df.countries.inc

# Joing cyclone affected data with country shape file and malaria incidence data
bi.map <- left_join(df.countries.inc, df.jing, by = join_by(ISO2 == ISO)) %>% 
  #Cleaning up data files
  mutate(avg_person_day = ifelse(is.na(avg_person_day), 0, avg_person_day)) %>%
  #Dropping Antartica
  filter(ISO2 != "AQ") %>%
  #creating my own bins for bivariate plotting
  #Breaks: Malaria: 0 (non-endemic), 0 to 10 (low), 10 to 100 (moderate), 100 to max (high)
  mutate(malaria.bin = cut(inc,            breaks = c(0, 0.0001,   10,    100, max(inc)),            include.lowest = TRUE)) %>%
  #Breaks: Storms: 0 (no storms), 0 to 80,0000 (low), 80,000 to 520,000 (moderate), 520,000 to max (high)
  mutate(storm.bin   = cut(avg_person_day, breaks = c(0, 0.001, 80000, 520000, max(avg_person_day)), include.lowest = TRUE))

#### Making bivariate map
#Calling bi_class to create a column with codes for bins for plotting
bi.map.c <- bi_class(bi.map, x = malaria.bin, y = storm.bin, style = "quantile", dim = 4)
sf.bi.map.c <- st_as_sf(x = bi.map.c, sf_column_name = "geometry") %>% st_set_crs(4326)

#Creating the map
p.bi.map <- ggplot() +
  #geom_sf(data = countries) +
  #Plotting, filling by bi_class. Dropping legend to make a better legend below
  geom_sf(data = sf.bi.map.c, mapping = aes(fill = bi_class), color = "white", size = 0.1, show.legend = FALSE) +
  bi_scale_fill(pal = "DkViolet2", dim = 4) +
  bi_theme()
#Making the legend
p.bi.map.legend <- bi_legend(pal = "DkViolet2",
                             dim = 4,
                             ylab = "Higher Cyclone Exposure",
                             xlab = "Higher Malaria Incidence",
                             size = 8)
#Plotting map + legend
p.bi.map + p.bi.map.legend + plot_layout(nrow = 1, widths = c(7, 2))



########################################################################################################################################
## Figure 1C ## Scatter plot of malaria incidence vs cyclone activity at national level
########################################################################################################################################

# Using ibtracs data

#Converting storm tracks line list data (dfa) to sf
sf.storm.points <- st_as_sf(dfa, coords = c("LON", "LAT"))
st_crs(sf.storm.points) <- 4326

#Using st_intersection to find storms track points that pass within buffer of malaria endemic countries
f.storm.point.container <- function(v.country_names, sf.points, sf.country_polygons){
  l.sf <- list(NA)
  for(i in 1:length(v.country_names)){
    contained_points_i <- st_contains(sf.country_polygons %>% filter(COUNTRY == v.country_names[i]), sf.points)
    sf_i <- sf.points[contained_points_i[[1]], ]
    sf_i <- sf_i %>% mutate(COUNTRY = v.country_names[i])
    l.sf[[i]] <- sf_i
  }
  df.sf <- bind_rows(l.sf)
  return(df.sf)
}

#note: takes approx. 5 minutes to run
sf.points.END <- f.storm.point.container(v.country_names = sf.END.buff$COUNTRY, 
                                         sf.points = sf.storm.points, 
                                         sf.country_polygons = sf.END.buff)

#Test plotting
ggplot() +
  #Base layer: gobal map from country shape files, dropping Antartica
  geom_sf(data = countries %>% filter(ISO != "AQ"),
          fill = "grey80", alpha = 0.2, color = "grey80", linewidth = 0.1, show.legend = FALSE) +
  #endemic hits
  geom_sf(data = sf.points.END, aes(color = COUNTRY), alpha = 0.1) + 
  coord_sf() +
  theme_bw() +
  theme(legend.position = "none")

#Tidying the data
sf.hits.END <- sf.points.END %>% 
  #dropping the n column for clarity (will be counting data points on subsequent groups)
  dplyr::select(-n_points) %>%
  #convert date to date format with lubridate
  mutate(date = ymd(substr(ISO_TIME, 1, 10))) %>%
  #Finding the max wind speed on storm approach (within 60nm of land)
  group_by(COUNTRY, name.yr) %>% 
  mutate(n = n(), n_NAs = sum(is.na(wind.speed)), n_wind.points = n()-n_NAs) %>%
  #For storms with wind speed data, find point with max wind
  mutate(max.wind.speed = ifelse(n_wind.points > 0, max(wind.speed, na.rm = TRUE), NA)) %>% ungroup() %>%
  filter((wind.speed == max.wind.speed) %>% replace_na(TRUE)) %>%
  #When multiple points have the same max wind, choose the one closest to land and earliest date
  group_by(COUNTRY, name.yr) %>%
  mutate(min.DIST2LAND = min(DIST2LAND)) %>% ungroup() %>%
  filter(DIST2LAND == min.DIST2LAND) %>%
  group_by(COUNTRY, name.yr) %>%
  mutate(min.date = min(date)) %>% ungroup() %>%
  filter(date == min.date) %>%
  #When max wind was the same at multiple points on earliest date, arbitrarily keep earliest data point
  group_by(COUNTRY, name.yr) %>%
  arrange(ISO_TIME) %>%
  filter(row_number()==1) %>%
  mutate(n.points = n()) %>% ungroup() %>%
  #Adding the year of the hit
  mutate(year = year(date)) %>%
  #Adding storm severity category
  mutate(storm.cat = case_when(
    max.wind.speed < 62                              ~ "TD",
    max.wind.speed > 62    & max.wind.speed <= 118.5 ~ "TS",
    max.wind.speed > 118.5 & max.wind.speed <= 153   ~ "H1",
    max.wind.speed > 153   & max.wind.speed <= 177   ~ "H2",
    max.wind.speed > 177   & max.wind.speed <= 208   ~ "H3",
    max.wind.speed > 208   & max.wind.speed <= 251   ~ "H4",
    max.wind.speed > 251                             ~ "H5",
    is.na(max.wind.speed) ~ "ND")) %>%
  #setting factor levels
  mutate(storm.cat = factor(storm.cat, levels = c("ND", "TD", "TS", "H1", "H2", "H3", "H4", "H5"))) %>%
  #arranging by country
  arrange(COUNTRY)

#### Tidying and exporting for ease of use
df.out.sf.hits.END <- sf.hits.END %>% st_drop_geometry() %>% 
  dplyr::select(COUNTRY, SEASON:NAME, name.yr, date, year, n_wind.points, max.wind.speed, storm.cat)
# write_csv(df.out.sf.hits.END, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_ibtracs_cyclones/processed/endemic_country_cyclone_hit_points.csv")



#########################################################################################################
# Extracting summary stats

#total number of storm hits among malaria endemic countries
length(df.out.sf.hits.END$COUNTRY)
#total number of storms that hit at least one malaria endemic country
length(unique(df.out.sf.hits.END$name.yr))
#total for storms with wind speed data and above Cat 1 hurricane level winds
length(df.out.sf.hits.END$COUNTRY[df.out.sf.hits.END$storm.cat %in% c("H1", "H2", "H3", "H4", "H5")])
#total for 'major' storms above Cat 3
length(df.out.sf.hits.END$COUNTRY[df.out.sf.hits.END$storm.cat %in% c("H3", "H4", "H5")])

#summarizing by country
df.all.storms <- sf.hits.END %>% st_drop_geometry() %>% group_by(COUNTRY) %>% summarize(n.ib.storms = n())
df.h1.plus    <- sf.hits.END %>% st_drop_geometry() %>% filter(storm.cat %in% c("H1", "H2", "H3", "H4", "H5")) %>% group_by(COUNTRY) %>% summarize(n.ib.storms = n())

#joining with malaria incidence data and malaria case count totals
df.storm.mal.all    <- df.all.storms %>% left_join(df.countries.inc, by = join_by(COUNTRY))
df.storm.mal.h1plus <- df.h1.plus    %>% left_join(df.countries.inc, by = join_by(COUNTRY))



#########################################################################################################
#Plotting

#Option 1
#Plotting malaria incidence vs storm freq (from ibtracs)
#All storms (tropical depression, etc)
#Size: Number of storms
#Color: Malaria incidence
p.m.c.scatter.ibtracs <- df.storm.mal.all %>% 
  #Labeling by 3 letter country code
  ggplot(aes(x = n.ib.storms, y = inc, 
             color = inc, 
             label = ISO3)) + 
  geom_point(aes(size=n.ib.storms)) +
  #Adding labels, adjusting padding to get labels to be clear
  geom_text_repel(#data = df2 %>% filter(num_storms > 0), 
    aes(color=inc),
    box.padding   = 1.8,
    max.overlaps  = 50,
    #making lines connecting dots to labels thinner
    segment.size  = 0.25) +
  #Manually adjusting scale, using log10 scale
  scale_y_continuous(limits = c(0.001, 400),
                     trans  = 'log10', 
                     breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                     labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
  scale_x_continuous(limits = c(1, 450),
                     trans  = 'log10', 
                     breaks = c(1, 10, 100),
                     labels = c(1, 10, 100)) +
  scale_color_viridis_c(option = "rocket",
                        begin = 0,
                        end = 0.7,
                        trans = "log2",
                        breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                        labels = label_number(accuracy = 0.1),
                        name = "Malaria\nincidence\n(per 1,0000 popn)") +
  #Adjusting the range in point size to highlight difference
  scale_size_continuous(name = "Storm count",
                        range = c(1,18)) +
  xlab("Number of storms since 1980 in IBTrACS data set") + ylab("Malaria incidence (per 1,000 population)") +
  #THEME: Making axes larger for legibility and removing some grid lines
  theme_bw() +
  theme(axis.text.x = element_text(size=13, family = "Arial"),
        axis.text.y = element_text(size=13, family = "Arial"),
        axis.title  = element_text(size=16),
        panel.grid.minor = element_blank())
p.m.c.scatter.ibtracs


#Option 2:
#Size: Number of cases in 2021
p.m.c.scatter.cases <- df.storm.mal.all %>% 
  ggplot(aes(x = n.ib.storms, y = inc, 
             label = ISO3)) + 
  geom_point(aes(size = cases), color = "#D2293D") +
  geom_text_repel(
    color = "#D2293D",
    box.padding   = 1.5,
    max.overlaps  = 50,
    segment.size  = 0.25) +
  #Manually adjusting scale, using log10 scale
  scale_y_continuous(limits = c(0.001, 400),
                     trans  = 'log10', 
                     breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                     labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
  scale_x_continuous(limits = c(1, 450),
                     trans  = 'log10', 
                     breaks = c(1, 10, 100),
                     labels = c(1, 10, 100)) +
  # scale_color_viridis_c(option = "rocket",
  #                       begin = 0,
  #                       end = 0.7,
  #                       trans = "log2",
  #                       breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
  #                       labels = label_number(accuracy = 0.1),
  #                       name = "Malaria\nincidence\n(per 1,0000 popn)") +
  #Adjusting the range in point size to highlight difference
  scale_size_continuous(name = "Estimated number\nof malaria cases nationally",
                        range = c(1,22),
                        breaks = c(100, 100000, 1000000, 5000000, 10000000),
                        #trans = "log10",
                        labels = label_comma(accuracy = 1)) +
  xlab("Number of storms since 1980 in IBTrACS data set") + ylab("Malaria incidence (per 1,000 population)") +
  #THEME: Making axes larger for legibility and removing some grid lines
  theme_bw() +
  theme(axis.text.x = element_text(size=13, family = "Arial"),
        axis.text.y = element_text(size=13, family = "Arial"),
        axis.title  = element_text(size=16),
        panel.grid.minor = element_blank())
p.m.c.scatter.cases


#########################################################################################################
#Extracting summary stats for case counts

#Case counts among cyclone vulnerable malaria endemic countries
df.case_counts <- df.storm.mal.all %>%
  filter(inc > 0) %>%
  filter(n.ib.storms > 0)
  #filter(n.ib.storms > 14) 
sum(df.case_counts$cases)

#### Tidying and exporting for ease of use
df.case_counts.out <- df.storm.mal.all %>% st_drop_geometry() %>% 
  dplyr::select(COUNTRY, ISO3, n.ib.storms, inc, cases) %>%
  arrange(desc(n.ib.storms), desc(cases))
sum(df.case_counts.out$n.ib.storms)
# write_csv(df.case_counts.out, "/Users/blrice/Library/CloudStorage/Dropbox/R DROPBOX/2022_MNJ_TAZO/data/data_ibtracs_cyclones/processed/cyclone_and_case_counts.csv")

#Plotting Madagascar
p.mada.cycs_per_yr <- df.out.sf.hits.END %>% filter(COUNTRY == "Madagascar") %>%
  group_by(year, storm.cat) %>%
  summarize(n = n()) %>% ungroup() %>%
  #can drop 2024 as data is not yet complete
  #filter(year < 2024) %>%
  ggplot(aes(x = year, y = n, fill = storm.cat)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_viridis_d(option = "magma",
                       name = "Storm\nCategory",
                       direction = 1,
                       end = 0.90) +
  xlab("Year") +
  ylab("Number of storms per year for Madagascar") +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", angle = 0, hjust = 0.5, vjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(family = "Arial", size = 14),
        strip.background = element_rect(fill = "white"))
p.mada.cycs_per_yr















########################################################################################################################################
## Supplemental Figure S1B-C ## Cyclone tracks in Africa
########################################################################################################################################


#List of African countries
african_countries <- c("Angola", "Burundi", "Benin", "Burkina Faso", "Botswana", "Central African Republic", "Côte d'Ivoire", "Cameroon", "Congo DRC", "Congo", "Comoros", "Cabo Verde", "Djibouti", "Algeria", "Egypt", "Eritrea", "Ethiopia", "Réunion", "Mayotte", "Gabon", "Ghana", "Guinea", "Gambia", "Guinea-Bissau", "Equatorial Guinea", "Kenya", "Liberia", "Libya", "Lesotho", "South Africa", 
                       "Morocco", "Madagascar", "Mali", "Mozambique", "Mauritania", "Mauritius", "Malawi", "Namibia", "Niger", "Nigeria", "Rwanda", "Sudan", "South Sudan", "Senegal", "Sierra Leone", "Somalia", "Sao Tome and Principe", "Eswatini", "Seychelles", "Chad", "Togo", "Tunisia", "Tanzania", "Uganda", "Zambia", "Zimbabwe")
#Dropping Antartica
df.AF <- bi.map %>% filter(ISO2 != "AQ") %>% 
  #Adding a flag for Africa and filtering
  mutate(africa_flag = ifelse(COUNTRY %in% african_countries, 1, 0)) %>%
  filter(africa_flag == 1)

sf.AF <- st_as_sf(x = df.AF, sf_column_name = "geometry") %>% st_set_crs(4326)


#Making a 60nm = 111km buffer around African countries
sf.AF.buff <- st_buffer(sf.AF, dist = 111*1000)

#Using st_intersection to find storms tracks that pass within buffer of African countries
int.sf.AF = st_intersection(sfa.t, sf.AF.buff)
#Subsetting storm track data to those that pass within buffer
sfa.AF <- sfa.t %>% filter(name.yr %in% int.sf.AF$name.yr) 
dfa.s <- dfa %>% group_by(name.yr) %>% summarize(SEASON = mean(SEASON)) %>% filter(name.yr %in% int.sf.AF$name.yr) %>% ungroup()
sfa.AF <- sfa.AF %>% full_join(dfa.s, by = join_by(name.yr)) 

p.AF <- ggplot() + 
  #Base layer: African countries
  geom_sf(data = sf.AF,
          fill = "black", alpha = 0.85, color = "white", linewidth = 0) +
  #Countries with 60 nm buffer
  #geom_sf(data = sf.AF.buff %>% filter(ISO != "AQ"), alpha = 0, color = "grey30") +
  #Cyclone tracks
  geom_sf(data = sfa.AF, 
          aes(color = SEASON),
          alpha = 0.5, lwd = 2) +
  #Adding a layer of country borders so Madagascar is visible
  geom_sf(data = sf.AF,
          fill = "black", alpha = 0.3, color = "black", linewidth = 0.5) +
  scale_colour_viridis_c(option = "mako",
                         breaks = c(2024, 2014, 2004, 1994, 1984),
                         name = "Cyclone Year") +
  #Set coords
  coord_sf() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))
p.AF


#### Cyclone impacts by African country

df.AF.hits <- tibble(name.yr = int.sf.AF$name.yr,
                     COUNTRY = int.sf.AF$COUNTRY) %>%
  group_by(COUNTRY) %>% summarize(n = n()) %>% 
  mutate(prop = n/sum(n)) %>%
  arrange(n, COUNTRY) 
df.AF.hits <- df.AF.hits %>% 
  mutate(COUNTRY = factor(COUNTRY, levels = df.AF.hits$COUNTRY))

p.AF.hits <- df.AF.hits %>% 
  ggplot(aes(x = COUNTRY, y = n)) +
  geom_bar(stat = "identity", fill = "#33608C", alpha = 0.6) +
  xlab(NULL) + ylab("Number of tropical storm impacts since 1980") +
  theme_bw() +
  theme(axis.text.x = element_text(family = "Arial", angle = 90, hjust = 1, vjust = 0.5, size = 13),
        panel.grid = element_blank(),
        axis.title = element_text(family = "Arial", size = 16))
p.AF.hits


df.AF.mean.hits <- tibble(name.yr = int.sf.AF$name.yr,
                          COUNTRY = int.sf.AF$COUNTRY,
                          year = as.numeric(substr(name.yr, 1, 4))) %>%
  group_by(COUNTRY, year) %>% summarize(n = n(), .groups = "drop") %>%
  complete(COUNTRY, year, fill = list(n = 0)) %>%
  group_by(COUNTRY) %>% mutate(cum.n = cumsum(n)) %>%
  mutate(max.count = max(cum.n)) %>% ungroup() %>%
  arrange(max.count)

df.AF.mean.hits <- df.AF.mean.hits %>% 
  mutate(COUNTRY = factor(COUNTRY, levels = unique(df.AF.mean.hits$COUNTRY)))

p.AF.hits.cumulative <- df.AF.mean.hits %>% 
  ggplot(aes(x = year, y = cum.n, group = COUNTRY)) +
  geom_path(color = "#33608C", alpha = 0.6, linewidth = 0.9) +
  geom_text_repel(data = df.AF.mean.hits %>% filter(year == 2024) %>% filter(max.count > 10), 
                  aes(x = year, y = cum.n, label = COUNTRY),
                  color = "#33608C", 
                  min.segment.length = 0,
                  nudge_x = 5, 
                  hjust = 0,
                  segment.size = 0.2,
                  direction = "y",
                  point.padding = 0.3,
                  segment.curvature = -1e-20,
                  segment.ncp = 3,
                  #segment.angle = 90
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)), breaks = c(1980, 1990, 2000, 2010, 2020)) +
  xlab("Year") + ylab("Cumulative number of recorded tropical cyclones since 1980") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(family = "Arial", size = 12),
        axis.title = element_text(family = "Arial", size = 16))

p.AF.hits + p.AF.hits.cumulative










########################################################################################################################################
## 1d ## Cyclones in SE Madagascar and Mananjary district
########################################################################################################################################

#Making a 60nm = 111km buffer around MNJ
sf.MNJ.buff <- st_buffer(mg.adm2 %>% filter(ADM2_EN == "Mananjary"), dist = 111*1000)

#Using st_intersection to find storms tracks that pass within buffer of MNJ
int.sf.MNJ = st_intersection(sfa.t, sf.MNJ.buff)
#Subsetting storm track data to those that pass within buffer
sfa.MNJ <- sfa.t %>% filter(name.yr %in% int.sf.MNJ$name.yr) 
dfa.s.MNJ <- dfa %>% group_by(name.yr) %>% summarize(SEASON = mean(SEASON)) %>% filter(name.yr %in% int.sf.MNJ$name.yr) %>% ungroup()
sfa.MNJ <- sfa.MNJ %>% full_join(dfa.s.MNJ, by = join_by(name.yr)) 


p.MNJ <- ggplot() + 
  #Base layer: Madagascar
  geom_sf(data = mg.adm0, fill = "black", alpha = 0.9, color = "black") +
  
  #All Africa contacting Cyclone tracks
  geom_sf(data = sfa.AF %>% filter(! name.yr %in% sfa.MNJ$name.yr), color = "grey80", alpha = 0.15, linewidth = 3) +
  
  #MNJ contacting Cyclone tracks
  geom_sf(data = sfa.MNJ %>% filter(name.yr != "2022_2022025S11091_BATSIRAI") %>% filter(name.yr != "2023_2023036S12118_FREDDY") %>% filter(name.yr != "2023_2023061S22036_FREDDY"), 
          color = "#33608C", linewidth = 3, alpha = 0.5) +
  
  #Cyclone Batsirai
  ##Points
  geom_point(data = dfa %>% filter(name.yr=="2022_2022025S11091_BATSIRAI"), aes(x=LON, y=LAT, group=name.yr, color=wind.speed), size=5.3, alpha=0.9) +
  ##Wide path with wind speed
  geom_path(data = dfa  %>% filter(name.yr=="2022_2022025S11091_BATSIRAI"), aes(x=LON, y=LAT, group=name.yr, color=wind.speed), linewidth=6, alpha=0.9) +
  ##Line
  geom_path(data = dfa  %>% filter(name.yr=="2022_2022025S11091_BATSIRAI"), aes(x=LON, y=LAT, group=name.yr), linewidth=0.5, alpha=0.9, color="grey30") +
  
  #Cyclone Freddy part 1
  ##Points
  geom_point(data = dfa %>% filter(name.yr=="2023_2023036S12118_FREDDY"), aes(x=LON, y=LAT, group=name.yr, color=wind.speed), size=5.3, alpha=0.9) +
  ##Wide path with wind speed
  geom_path(data = dfa  %>% filter(name.yr=="2023_2023036S12118_FREDDY"), aes(x=LON, y=LAT, group=name.yr, color=wind.speed), linewidth=6, alpha=0.9) +
  ##Line
  geom_path(data = dfa  %>% filter(name.yr=="2023_2023036S12118_FREDDY"), aes(x=LON, y=LAT, group=name.yr), linewidth=0.5,alpha=0.9, color="grey30") +
  
  #Cyclone Freddy part 2
  ##Points
  #geom_point(data = dfa %>% filter(name.yr=="2023_2023061S22036_FREDDY"), aes(x=LON, y=LAT, group=name.yr, color=wind.speed), size=5.3, alpha=0.9) +
  ##Wide path with wind speed
  #geom_path(data =  dfa %>% filter(name.yr=="2023_2023061S22036_FREDDY"), aes(x=LON, y=LAT, group=name.yr, color=wind.speed), linewidth=6, alpha=0.9) +
  ##Line
  #geom_path(data =  dfa %>% filter(name.yr=="2023_2023061S22036_FREDDY"), aes(x=LON, y=LAT, group=name.yr), linewidth=0.5, alpha=0.9, color="grey30") +
  
  #MNJ District
  geom_sf(data = mg.adm2 %>% filter(ADM2_EN == "Mananjary"),
          fill = "black", alpha = 0, color = "white", linewidth = 0.4) +
  scale_color_viridis_c(option = "magma",
                        name = "Wind Speed\n(km/h)",
                        end = 0.9,
                        limits = c(0, 260)) +
  #Set coords
  coord_sf(xlim = c(38.5, 75), ylim = c(-11, -27)) +
  #coord_sf(xlim = c(25, 125), ylim = c(-11, -32)) +
  
  theme(panel.grid = element_blank(),
        legend.position = "right", legend.justification = "bottom",
        #legend.position = "none",
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))
p.MNJ









