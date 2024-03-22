
library("sf")

sdn_map <- readRDS(url("https://github.com/jalilian/CEASE/raw/main/Sudan/sdn_map.rds"))

source("https://raw.githubusercontent.com/jalilian/CEASE/main/Ethiopia/codes/get_land_covars_africa.R")

sdn_covars <- get_covars(sdn_map, path="~/Downloads/Africa_covars/")
