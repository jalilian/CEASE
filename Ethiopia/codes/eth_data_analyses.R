
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

# merge Woreds in each zone
eth_map <- eth_map %>%
  group_by(ADM1_EN, ADM2_EN) %>%
  summarise(geometry=st_union(geometry),
            Total_pop=sum(Total)) %>%
  ungroup() %>%
  rename(RegionName=ADM1_EN, ZoneName=ADM2_EN)
# =========================================================
# missing data patterns

# number of zones with recorded cases by date
tmp1 <- data.frame(Date=seq(as.Date("2013-01-07"), 
                            as.Date("2022-08-15"), 
                            by=7)) %>%
  left_join(eth_data %>% count(Date),
            by = join_by(Date)) %>%
  mutate(n=if_else(is.na(n), 0, n))

tmp1 %>%
  ggplot(aes(x=Date, y=(1 - n / 92))) +
  geom_line() +
  labs(y="percentage of missing weekly malaria records by date") +
  scale_y_continuous(labels=scales::percent)  +
  theme_light()

library("tseries")
# the augmented Dickey-Fuller (ADF) test for stationarity (unit root)
adf.test(tmp1 %>% pull(n))
# Kwiatkowski-Phillips-Schmidt-Shin (KPSS) for stationarity (constant mean)
kpss.test(tmp1 %>% pull(n))

# number of weeks with recorded cases by zones
tmp2 <- eth_map %>%
  left_join(eth_data %>% 
              count(ZoneName, Year) %>%
              pivot_wider(names_from=Year, values_from=n),
            by=join_by(ZoneName)) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols=`2013`:`2022`, 
               names_to="Year", 
               values_to ="n") 

tmp2 %>%
  ggplot(aes(fill=1 - n/374)) + 
  geom_sf() +
  scale_fill_distiller(name="missing", 
                       labels=scales::percent,
                       palette = "Reds", direction=1)+
  facet_wrap(~Year, ncol=4) +
  theme_light() +
  theme(legend.position='bottom')

library("spdep")

for (year in c(2013:2019, 2022))
{
  W <- nb2listw(poly2nb(tmp2 %>% filter(Year == year)), 
                style = "B")
  # Moran's I test for spatial autocorrelation
  print(moran.test(tmp2 %>% filter(Year == year) %>% pull(n), 
             listw = W))
}

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
  left_join(eth_data %>% 
              group_by(ZoneName, Year) %>%
              summarise(Total_confirmed=
                          sum(`Total Malaria Confirmed and Clinical`)) %>%
              pivot_wider(names_from=Year, 
                        values_from=Total_confirmed),
            by=join_by(ZoneName)) %>% 
  pivot_longer(cols=`2013`:`2022`, 
               names_to="Year", 
               values_to ="Total_confirmed") %>%
  # account for population changes by considering constant 2.67% population growth rate
  mutate(Total_pop2=
           Total_pop * ( (1 - 2.67 / 100)^(2022 - as.numeric(Year)) )) %>%
  mutate(rate=Total_confirmed / Total_pop2 * 1e5) %>%
  ggplot(aes(fill=rate)) + 
  geom_sf() +
  scale_fill_distiller(name="confirmed", palette = "Reds", direction=1)+
  facet_wrap(~Year, ncol=4) +
  theme_light() +
  theme(legend.position='bottom',
        legend.key.width = unit(3, 'cm'))


