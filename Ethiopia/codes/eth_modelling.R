
library("tidyverse")

# =========================================================
# path to the malaria data file
data_path <- "~/Downloads/Ethiopia/"

# reading the malaria data
eth_data <- readRDS(paste0(data_path, "eth_data.rds"))

# =========================================================

# create spatial and temporal indices
eth_data <- eth_data %>%
  left_join(
    expand_grid(Year=2013:2022, 
                Epidemic_Week=1:52) %>%
      mutate(idx_time=1:520),
    by = join_by(Year, Epidemic_Week)
  ) %>%
  left_join(
    eth_data %>% 
      count(ZoneName2) %>% 
      select(-n) %>% 
      mutate(idx_zone=1:90),
    by = join_by(ZoneName2)
  )


