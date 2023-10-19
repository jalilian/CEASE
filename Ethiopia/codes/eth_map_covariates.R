
library("sf")
library("tidyverse")
library("terra")

# =========================================================
# read map data
eth_map <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))

# =========================================================
# function to read and extract raster data
extract_raster_data <- function(downloaded_file,
           data_file,
           temp_dir="/tmp/eth_raster/")
{
    # create a temporary directory to extract the downloaded file
    dir.create(temp_dir)
    # extract the downloaded file
    unzip(zipfile=downloaded_file, exdir=temp_dir)
    
    # read the raster data using the 'terra' package
    rast(paste0(temp_dir, data_file))
}

# =========================================================
# read the Gridded Population of the World (GPW), v4
# provided by Socioeconomic Data and Applications Center, NASA
# https://sedac.ciesin.columbia.edu/data/collection/gpw-v4
pop_density <- 
  extract_raster_data(
    "~/Downloads/gpw-v4-population-density-rev11_2020_30_sec_tif.zip",
    "gpw_v4_population_density_rev11_2020_30_sec.tif"
  )

# read S2 Prototype Land Cover 20m Map of Africa 2016
# provided by CCI Land Cover (LC) team, European Space Agency
# https://2016africalandcover20m.esrin.esa.int/
africa_land_cover <- 
  extract_raster_data(
    "~/Downloads/ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.zip",
    "ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif"
  )

# =========================================================
# plot population density and land cover for the whole country
library("tidyterra")
ggplot() +
  geom_spatraster(data=mask(crop(pop_density, eth_map),
                            eth_map),
                  aes(fill=gpw_v4_population_density_rev11_2020_30_sec))
ggplot() +
  geom_spatraster(data=mask(crop(africa_land_cover, eth_map),
                            eth_map),
                  aes(fill=`ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0`))

# =========================================================
# extract population density and land cover data for map units

# merge Woreds in each zone
map <- 
  eth_map %>%
  group_by(ADM2_EN) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# extract covariates
covar_fun <- function(map_unit)
{
  cp1 <- terra::crop(pop_density, map_unit)
  cp1 <- terra::mask(cp1, map_unit)
  cp1 <- values(cp1)
  cp1 <- cp1[!is.na(cp1)]
  
  cp2 <- terra::crop(africa_land_cover, map_unit)
  cp2 <- terra::mask(cp2, map_unit)
  cp2 <- values(cp2)
  cp2 <- cp2[!is.na(cp2)]
  cp2 <- sort(table(cp2), decreasing=TRUE)
  cp2 <- round(100 * cp2[1] / sum(cp2), 2)
  
  return(c(popdens_mean=mean(cp1), 
           popdens_sd=sd(cp1), 
           landcov_mode=as.integer(names(cp2)),
           landcov_percent=unname(cp2)))
}

covars <- 
  parallel::mclapply(split(map, 1:nrow(map)),
                     covar_fun, 
                     mc.cores=1)

covars <- 
  bind_cols(ADM2_EN=
              unique(eth_map %>% 
                       pull(ADM2_EN)), 
            bind_rows(covars))
