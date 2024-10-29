
library("sf")
library("tidyverse")

# =========================================================

# read the pre-proceed malaria data
eth_data <- 
  readRDS("~/Downloads/Ethiopia/eth_data.rds")

# remove all records of week 53 of the year 2015
eth_data %>%
  filter(is.na(Date))

# data from a shapefile provided by OCHA
eth_map <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))

# merge Woreds in each zone
eth_map <- eth_map %>%
  mutate(ADM2_EN=
           case_match(ADM2_EN,
                      "East Bale" ~ "Bale",
                      "Dire Dawa rural" ~ "Dire Dawa urban",
                      .default=ADM2_EN)) %>%
  group_by(ADM1_EN, ADM2_EN) %>%
  summarise(geometry=st_union(geometry),
            Total_pop=sum(Total)) %>%
  ungroup() %>%
  rename(RegionName=ADM1_EN, ZoneName=ADM2_EN)
# =========================================================
# missing data patterns

tmp <- expand_grid(
  ZoneName=eth_map %>% pull(ZoneName),
  eth_data %>% 
    count(Year, Epidemic_Week) %>% select(-n)
  )

# number of missing zones by date
tmp <- tmp %>%
  left_join(
    eth_data %>% 
      select(ZoneName, Year, Epidemic_Week) %>%
      mutate(missing=0),
    by = c("ZoneName", "Year", "Epidemic_Week")
  ) %>%
  mutate(missing = if_else(is.na(missing), 1, 0))

# by year
tmp %>%
  group_by(Year) %>%
  summarise(n=100 * mean(missing))

# by Zone
tmp %>%
  group_by(ZoneName) %>%
  summarise(n=100 * mean(missing))

tmp1 <- tmp %>% 
  mutate(Date=ymd(paste(Year, "01", "01", sep="-")) + 
           weeks(Epidemic_Week)) %>%
  group_by(Date) %>% 
  summarise(n=100 * mean(missing))

g1 <- tmp1 %>%
  ggplot(aes(x=Date, y=n)) + 
  geom_line() +
  stat_smooth(se=FALSE) +
  labs(x="date", y="percentage of weekly missing zones") +
  theme_light()

library("tseries")
# the augmented Dickey-Fuller (ADF) test for stationarity (unit root)
adf.test(tmp1 %>% pull(n))
# Kwiatkowski-Phillips-Schmidt-Shin (KPSS) for stationarity (constant mean)
kpss.test(tmp1 %>% pull(n))

# number of missing weeks by zone
tmp2 <- tmp %>% 
  group_by(ZoneName, Year) %>% 
  summarise(n=100 * mean(missing)) %>%
  ungroup()


g2 <- eth_map %>% 
  left_join(tmp2 %>% 
              group_by(ZoneName) %>%
              summarise(n=mean(n)), 
            by=join_by(ZoneName)) %>%
  ggplot(aes(fill= n)) + 
  geom_sf() +
  scale_fill_distiller(name="missing", 
                       labels=scales::percent,
                       palette = "Reds", direction=1)+
  theme_light() +
  theme(legend.position='bottom')

library("ggpubr")
library("devEMF")
ggarrange(g1, g2, nrow=1)
ggsave(filename="Missing.emf", device=emf, width=15, height=7)

tmp2 <- eth_map %>% 
  left_join(tmp2, by=join_by(ZoneName))

tmp2 %>%
  ggplot(aes(fill= n)) + 
  geom_sf() +
  scale_fill_distiller(name="missing", 
                       labels=scales::percent,
                       palette = "Reds", direction=1)+
  facet_wrap(~Year, ncol=4) +
  theme_light() +
  theme(legend.position='bottom')


# Moran's I test for spatial autocorrelation
library("spdep")

# overall
W <- nb2listw(poly2nb(
  tmp2 %>% 
    group_by(ZoneName) %>% 
    summarise(miss=mean(miss))
))


moran.test(
  tmp2 %>% 
    group_by(ZoneName) %>% 
    summarise(miss=mean(miss)) %>% pull(miss),
  listw = W)

# by Year
for (year in c(2013:2019, 2022))
{
  W <- nb2listw(poly2nb(tmp2 %>% filter(Year == year)), 
                style = "B")
  # Moran's I test for spatial autocorrelation
  print(moran.test(tmp2 %>% filter(Year == year) %>% pull(miss), 
             listw = W))
}


# =========================================================
# confirmed cases 

# temporal pattern
tmp1 <- eth_data %>% 
  count(Year, Epidemic_Week) %>% select(-n) %>%
  left_join(eth_data %>% 
              group_by(Year, Epidemic_Week) %>% 
              summarise(Total_confirmed=sum(`Total Malaria Confirmed and Clinical`)),
            by = join_by(Year, Epidemic_Week)) %>% 
  mutate(Date=ymd(paste(Year, "01", "01", sep="-")) + 
           weeks(Epidemic_Week))
  
tmp1 <- data.frame(Date=seq(as.Date("2013-01-07"), 
                            as.Date("2022-08-15"), 
                            by=7)) %>%
  left_join(eth_data %>% 
              group_by(Date) %>% 
              summarise(Total_confirmed=sum(`Total Malaria Confirmed and Clinical`)),
            by = join_by(Date))
tmp1 %>%
  ggplot(aes(x=Date, y=Total_confirmed)) +
  geom_line() +
  stat_smooth() +
  labs(y="Total number of clinically confirmed malaria cases") +
  theme_light()

# check ACF and PACF 
library("forecast")
ggAcf(tmp1$Total_confirmed)
ggPacf(tmp1$Total_confirmed)
ggAcf(diff(tmp1$Total_confirmed))
ggPacf(diff(tmp1$Total_confirmed))
ggAcf(diff(tmp1$Total_confirmed, 2))
ggPacf(diff(tmp1$Total_confirmed, 2))

library("tseries")
adf.test(na.omit(tmp1$Total_confirmed[tmp1$Date < ymd("2019-12-31")]))
adf.test(na.omit(diff(tmp1$Total_confirmed[tmp1$Date < ymd("2019-12-31")])))
adf.test(na.omit(diff(tmp1$Total_confirmed[tmp1$Date < ymd("2019-12-31")], 2)))

# spatial pattern
eth_map %>%
  left_join(eth_data %>% 
              group_by(ZoneName, Year) %>%
              summarise(Total_confirmed=
                          sum(`Total Malaria Confirmed and Clinical`)) %>%
              pivot_wider(names_from=Year, 
                        values_from=Total_confirmed),
            by=join_by(ZoneName)) %>% 
  pivot_longer(cols=`2013`:`2023`, 
               names_to="Year", 
               values_to ="Total_confirmed") %>%
  # account for population changes by considering constant 2.67% population growth rate
  mutate(Total_pop2=
           Total_pop * ( (1 - 2.67 / 100)^(2022 - as.numeric(Year)) )) %>%
  mutate(rate=Total_confirmed / Total_pop2 * 1e4) %>%
  ggplot(aes(fill=rate)) + 
  geom_sf() +
  scale_fill_distiller(name="confirmed", palette = "Reds", direction=1)+
  facet_wrap(~Year, ncol=4) +
  theme_light() +
  theme(legend.position='bottom',
        legend.key.width = unit(3, 'cm'))


