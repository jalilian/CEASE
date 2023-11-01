
library("sf")
library("tidyverse")

# =========================================================
# path to the malaria data file
data_path <- "~/Downloads/Ethiopia/"

# reading the malaria data
eth_data <- readRDS(paste0(data_path, "eth_data.rds"))

# data from a shapefile provided by OCHA
eth_map <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))

# merge Woreds in each zone
eth_map <- eth_map %>%
  group_by(ADM1_EN, ADM2_EN) %>%
  summarise(geometry=st_union(geometry),
            Total=sum(Total)) %>%
  ungroup() %>%
  rename(RegionName=ADM1_EN, ZoneName=ADM2_EN)
# =========================================================

# remove all records of week 53 of the year 2015
eth_data <- eth_data %>%
  filter(!is.na(Date))

eth_data %>% 
  left_join(tibble(Date=seq(as.Date("2013-01-07"), 
                            as.Date("2022-08-15"), 
                            by=7)) %>%
              mutate(idx_time=1:length(Date)),
            by=join_by(Date)) %>%
  left_join(eth_map %>% 
              mutate(idx_zone=1:length(ZoneName)) %>%
              as_tibble() %>%
              select(RegionName, ZoneName, idx_zone),
            by=join_by(RegionName, ZoneName))



# =========================================================
library("INLA")