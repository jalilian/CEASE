
library("sf")
library("tidyverse")

# =========================================================

# read the pre-proceed malaria data
eth_data <- 
  readRDS("~/Downloads/Ethiopia/eth_data.rds")

# remove all records of week 53 of the year 2015
eth_data <- eth_data %>%
  filter(!is.na(Date))


eth_data %>%
  full_join(
    data.frame(Date=seq(as.Date("2013-01-07"), 
                        as.Date("2022-08-15"), 
                        by=7)),
    by = join_by(Date)
  )


eth_data %>% 
  count(Date) %>%
  ggplot(aes(x=Date, y=n)) +
  geom_line()
