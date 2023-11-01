
library("sf")
library("tidyverse")
library("terra")

# =========================================================
# read map data
eth_map <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))

# merge Woreds in each zone
eth_map <- eth_map %>%
  group_by(ADM2_EN) %>%
  summarise(geometry=st_union(geometry),
            Total=sum(Total)) %>%
  ungroup()

# =========================================================
# spatial covariates
# =========================================================
# function to read and extract raster data
extract_raster_covars <- function(downloaded_file,
                                  data_file,
                                  map,
                                  stat="mean",
                                  temp_dir="/tmp/eth_raster/")
{
  # create a temporary directory to extract the downloaded file
  if (!dir.exists(temp_dir))
  {
    dir.create(temp_dir)
  }
  # extract the downloaded file
  unzip(zipfile=downloaded_file, exdir=temp_dir)
  
  # read the raster data using the 'terra' package
  raster <- rast(paste0(temp_dir, data_file))
  
  out <- vector("list", length=nrow(map))
  for (i in 1:nrow(map))
  {
    cp <- terra::crop(raster, map[i, ], mask=TRUE, touches=TRUE)
    cp <- values(cp)
    cp <- cp[!is.na(cp)]
    
    out[[i]] <- switch(stat, mean={
      c(mean=mean(cp), sd=sd(cp))
    }, mode={
      cp <- sort(table(cp), decreasing=TRUE)
      cp <- round(100 * cp[1] / sum(cp), 2)
      c(mode=names(cp), percent=unname(cp))
    })
  }
  
  # clean up the temporary directory
  unlink(temp_dir, recursive=TRUE)
  
  return(bind_rows(out))
}

# read the Gridded Population of the World (GPW), v4
# resolution: 30 arc second
# provided by Socioeconomic Data and Applications Center, NASA
# https://sedac.ciesin.columbia.edu/data/collection/gpw-v4
pop_density <- 
  extract_raster_covars(
    paste0("~/Downloads/Africa_covars/",
           "gpw-v4-population-density-rev11_2020_30_sec_tif.zip"),
    "gpw_v4_population_density_rev11_2020_30_sec.tif",
    map=eth_map, stat="mean"
  )

# read S2 Prototype Land Cover 20m Map of Africa 2016
# resolution: 20 meter
# provided by CCI Land Cover (LC) team, European Space Agency
# https://2016africalandcover20m.esrin.esa.int/
land_cover <- 
  extract_raster_covars(
    paste0("~/Downloads/Africa_covars/",
    "ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.zip"),
    "ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif",
    map=eth_map, stat="mode"
  )

# Void-filled digital elevation mode of Africa, 2007 [ArcGRID]
# resolution: 15s resolution
# provided by World Wildlife Fund. Africa, U.S. Geological Survey
# https://maps.princeton.edu/catalog/stanford-jm998sr5227 
elevation <- 
  extract_raster_covars(
    paste0("~/Downloads/Africa_covars/",
    "data.zip"),
    "af_dem_15s/",
    map=eth_map, stat="mean"
  )

# combine raster covariates
spat_covars <-
  bind_cols(
    ADM2_EN=eth_map %>% pull(ADM2_EN), 
    pop_density %>% 
      setNames(paste0('pop_density_', names(.))),
    land_cover %>% 
      setNames(paste0('land_cover_', names(.))),
    elevation %>% 
      setNames(paste0('elevation_', names(.)))
  )

# save the extracted spatial covariates
saveRDS(spat_covars, file="spat_covars.rds")

# Continent of Africa: High Resolution Population Density Maps
# provided by Data for Good at Meta
# https://data.humdata.org/dataset/highresolutionpopulationdensitymaps


pop_dens_meta <- function(downloaded_file="popdens.zip",
                          temp_dir="/tmp/pop_dens_africa/")
{
  # create a temporary directory to extract the downloaded file
  if (!dir.exists(temp_dir))
  {
    dir.create(temp_dir)
  }
  download.file(
  "https://data.humdata.org/dataset/dbd7b22d-7426-4eb0-b3c4-faa29a87f44b/resource/7b3ef0ae-a37d-4a42-a2c9-6b111e592c2c/download/population_af_2018-10-01.zip", 
  destfile=paste0(temp_dir, downloaded_file)
  )
  # extract the downloaded file
  unzip(zipfile=paste0(temp_dir, downloaded_file), exdir=temp_dir)
  
  # read all raster files
  raster <- sapply(list.files(temp_dir, 
                              pattern=".tif$", 
                              full.names=TRUE), 
                   rast)

  raster <- sprc(raster)
  raster <- terra::crop(raster, ext(map))
  raster <- terra::mosaic(raster)
  
  out <- vector("list", length=nrow(map))
  for (i in 1:nrow(map))
  {
    cp <- terra::crop(raster, map[i, ], mask=TRUE, touches=TRUE)
    cp <- values(cp)
    cp <- cp[!is.na(cp)]
    
    out[[i]] <- switch(stat, mean={
      c(mean=mean(cp), sd=sd(cp))
    }, mode={
      cp <- sort(table(cp), decreasing=TRUE)
      cp <- round(100 * cp[1] / sum(cp), 2)
      c(mode=names(cp), percent=unname(cp))
    })
  }
  
  return(bind_rows(out))
}
