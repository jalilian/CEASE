
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

# save map data as an R object of class sf
saveRDS(eth_map, file="eth_map.rds")
