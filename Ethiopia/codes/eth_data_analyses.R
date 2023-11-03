
library("sf")
library("tidyverse")

# =========================================================

# read the pre-proceed malaria data
eth_data <- 
  readRDS("~/Downloads/Ethiopia/eth_data.rds")

# remove all records of week 53 of the year 2015
eth_data <- eth_data %>%
  filter(!is.na(Date))

# data from a shapefile provided by OCHA
eth_map <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))

# =========================================================
# missing data patterns

# number of zones with recorded cases by date
data.frame(Date=seq(as.Date("2013-01-07"), 
                    as.Date("2022-08-15"), 
                    by=7)) %>%
  left_join(eth_data %>% count(Date),
            by = join_by(Date)) %>%
  ggplot(aes(x=Date, y=(1 - n / 92))) +
  geom_line() +
  labs(y="percentage of missing weekly malaria records by date") +
  scale_y_continuous(labels=scales::percent)  +
  theme_light()


# number of weeks with recorded cases by zones
eth_map %>% 
  group_by(ADM2_EN) %>%
  summarise(geometry=st_union(geometry)) %>%
  ungroup() %>%
  rename(ZoneName=ADM2_EN) %>%
  left_join(eth_data %>% count(ZoneName),
            by=join_by(ZoneName)) %>%
  ggplot(aes(fill=1 - n/374)) + 
  geom_sf() +
  scale_fill_distiller(name="missing", 
                       labels=scales::percent,
                       palette = "Reds", direction=1)+
  theme_light()

# =========================================================
# confirmed cases 

# temporal pattern
data.frame(Date=seq(as.Date("2013-01-07"), 
                    as.Date("2022-08-15"), 
                    by=7)) %>%
  left_join(eth_data %>% 
              group_by(Date) %>% 
              summarise(Total_confirmed=sum(`Total Malaria Confirmed and Clinical`)),
            by = join_by(Date)) %>%
  ggplot(aes(x=Date, y=Total_confirmed)) +
  geom_line() +
  labs(y="Total number of clinically confirmed malaria cases") +
  theme_light()

# spatial pattern
eth_map %>% 
  group_by(ADM2_EN) %>%
  summarise(geometry=st_union(geometry)) %>%
  ungroup() %>%
  rename(ZoneName=ADM2_EN) %>%
  left_join(eth_data %>% 
              group_by(ZoneName) %>%
              summarise(Total_confirmed=sum(`Total Malaria Confirmed and Clinical`)),
            by=join_by(ZoneName)) %>%
  ggplot(aes(fill=Total_confirmed)) + 
  geom_sf() +
  scale_fill_distiller(name="confirmed", palette = "Reds", direction=1)+
  theme_light()

