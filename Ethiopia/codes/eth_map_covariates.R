
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



# =========================================================
library("ecmwfr")
library("ncdf4")
library("ncdf4.helpers")

# European Centre for Medium-Range Weather Forecasts (ECMWF) 
# Copernicus's Climate Data Store (CDS)
# https://cds.climate.copernicus.eu
# DOI: 10.24381/cds.e2161bac

# user credentials
user <- "****************"
cds.key <- "********************************"

# set secret ECMWF token
wf_set_key(user=user, key=cds.key, service="cds")

land_data_fun <- function(year, 
                          temp_dir="/tmp/eth_land/")
{
  # create a temporary directory to extract the downloaded file
  dir.create(temp_dir)
  setwd(temp_dir)
  
  # request for getting land data
  request <- list(
    # dataset name
    dataset_short_name = "reanalysis-era5-land",
    # climate variables 
    variable = c("10m_u_component_of_wind", 
                 "10m_v_component_of_wind", 
                 "leaf_area_index_high_vegetation", 
                 "leaf_area_index_low_vegetation", 
                 "skin_temperature", 
                 "surface_pressure",
                 "total_precipitation", 
                 "volumetric_soil_water_layer_1"),
    # temporal framework: year, month, day, hour
    year = as.character(year),
    month = sprintf("%02d", 1:12),
    day = sprintf("%02d", 1:31),
    time = "00:00",
    # geographical region
    #      North, West, South, East
    area = c(14.9, 32.9, 3.4, 48),
    # output file format
    format = "netcdf.zip",
    # output file name
    target = paste0("landvars_hourly_", year, ".zip")
  )
  
  # check the validity of a data request and login credentials
  wf_check_request(user=user, request=request)
  
  # download the data request
  wf_request(user=user, 
             request=request,
             transfer=TRUE, 
             path=getwd(),
             # waiting time for download to start
             time_out=3 * 60 * 60,
             verbose=TRUE)
  
  # extract downloaded Zip file
  unzip(zipfile=paste0(temp_dir, "landvars_hourly_", year, ".zip"), 
        exdir=paste0(temp_dir, year, "/"), 
        overwrite=TRUE)
  
  # open nc file containing the data
  nc_data <- nc_open(paste0(temp_dir, year, "/data.nc"))
  
  # extract longitude
  lon <- ncvar_get(nc_data, "longitude")
  # extract latitude
  lat <- ncvar_get(nc_data, "latitude")
  # extract date and time
  dt <- nc.get.time.series(nc_data)
  # list of names of data variables
  vars <- nc.get.variable.list(nc_data)
  
  # convert nc data to an R data.frame
  dat <- vector("list", length=length(vars))
  for (i in 1:length(vars))
  {
    vals <- ncvar_get(nc_data, vars[i])
    if (length(dim(vals)) > 3)
    {
      idx_3 <- 
        apply(vals, MARGIN=3, 
              function(x){ mean(is.na(x)) } 
        )
      idx_3 <- which.min(idx_3)
      vals <- vals[, , idx_3, ]
    }
    dimnames(vals) <- 
      list(longitude=1:length(lon), 
           latitude=1:length(lat),
           time=1:length(dt))
    
    dat[[i]] <- as.data.frame.table(vals)
    dat[[i]] <- dat[[i]] %>%
      mutate(longitude=lon[longitude],
             latitude=lat[latitude],
             time=dt[time]) %>%
      setNames(c(colnames(dat[[i]])[1:3], 
                 vars[i]))
  }
  
  dat %>% reduce(full_join, 
                 by = join_by(longitude, 
                              latitude, 
                              time))
}

land_covars <- 
  lapply(2013:2019, land_data_fun)

land_covars <- land_covars %>% 
  reduce(full_join)
