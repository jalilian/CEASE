
# Create a temporary directory to download map data
temp_dir <- "/tmp/eth"
dir.create(temp_dir)
setwd(temp_dir)

# download map data
map_link <- "https://data.humdata.org/dataset/cb58fa1f-687d-4cac-81a7-655ab1efb2d0/resource/63c4a9af-53a7-455b-a4d2-adcc22b48d28/download/eth_adm_csa_bofedb_2021_shp.zip"
download.file(url=map_link, destfile="eth_map.zip")

# extract Zip file
unzip(zipfile="eth_map.zip", exdir="eth_map/")

# read map data
library("sf")
eth_map <- read_sf("eth_map/",
                   layer="eth_admbnda_adm3_csa_bofedb_2021")

# plot map data
library("ggplot2")
ggplot(eth_map) + geom_sf()

# download population data
pop_link <- "https://data.humdata.org/dataset/3d9b037f-5112-4afd-92a7-190a9082bd80/resource/f82b20f1-8a76-46e9-ba9a-29e531f7af3c/download/eth_admpop_2023.xlsx"
download.file(url=pop_link, destfile="eth_pop.xlsx")

# read population data
library("readxl")
eth_pop <- read_excel("eth_pop.xlsx", 
                      sheet="ETH_admpop_adm3_2023",
                      range="A1:BD1085", 
                      col_names=TRUE, na="")

# change the names of columns representing administrative divisions
library("dplyr")
eth_pop <- eth_pop %>%
  rename(ADM3_EN=admin3Name_en, ADM3_PCODE=admin3Pcode,
         ADM2_EN=admin2Name_en, ADM2_PCODE=admin2Pcode,
         ADM1_EN=admin1Name_en, ADM1_PCODE=admin1Pcode,
         ADM0_EN=admin0Name_en, ADM0_PCODE=admin0Pcode)

# merge map and population data
eth_map <- eth_map %>% 
  left_join(eth_pop, by=c("ADM3_EN", "ADM3_PCODE",
                          "ADM2_EN", "ADM2_PCODE",
                          "ADM1_EN", "ADM1_PCODE",
                          "ADM0_EN", "ADM0_PCODE"))

# plot population data on the map
ggplot(eth_map) + geom_sf(aes(fill=Total))

# save map data as an R object of class sf
saveRDS(eth_map, file="eth_map.rds")
