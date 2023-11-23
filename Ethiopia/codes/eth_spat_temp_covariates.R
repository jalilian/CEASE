
library("tidyverse")
library("sf")
library("ecmwfr")
library("ncdf4")
library("ncdf4.helpers")

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
# spatiotemporal covariates
# =========================================================
# European Centre for Medium-Range Weather Forecasts (ECMWF) 
# Copernicus's Climate Data Store (CDS)
# https://cds.climate.copernicus.eu
# DOI: 10.24381/cds.e2161bac

# user credentials for ECMWF data access
user <- "****************"
cds.key <- "********************************"

# set secret ECMWF token
wf_set_key(user=user, key=cds.key, service="cds")

land_data_fun <- function(year, 
                          temp_dir="/tmp/eth_land/")
{
  # create a temporary directory to extract the downloaded file
  if (!dir.exists(temp_dir))
  {
    dir.create(temp_dir)
  }
  # set the working directory to the temporary directory
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
  
  # open netCDF file containing the data
  nc_data <- nc_open(paste0(temp_dir, year, "/data.nc"))
  
  # extract longitude
  lon <- ncvar_get(nc_data, "longitude")
  # extract latitude
  lat <- ncvar_get(nc_data, "latitude")
  # extract date and time
  dt <- nc.get.time.series(nc_data)
  # list of names of data variables
  vars <- nc.get.variable.list(nc_data)
  
  # create a spatial grid using longitude and latitude
  grid <- 
    expand_grid(lon_idx=1:length(lon),
                lat_idx=1:length(lat)) %>%
    rowwise() %>% 
    mutate(points=list(st_point(c(lon[lon_idx], 
                                  lat[lat_idx])))) %>%
    st_as_sf()
  
  # set the coordinate reference system for the grid
  st_crs(grid) <- st_crs(eth_map)
  
  # conduct a spatial join to determine points inside specific map regions
  grid <- 
    st_join(grid, eth_map, join=st_within) %>%
    filter(!is.na(ADM2_EN))

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
    
    # define dimension names and indices
    dimnames(vals) <- 
      list(lon_idx=1:length(lon), 
           lat_idx=1:length(lat),
           time_idx=1:length(dt))
    
    # convert data array to a data frame
    dat[[i]] <- as.data.frame.table(vals) %>%
      mutate(lon_idx=as.integer(lon_idx),
             lat_idx=as.integer(lat_idx),
             time_idx=as.integer(time_idx))
    
    # aggregate by map regions
    # determine grid points inside map regions
    dat[[i]] <- dat[[i]] %>%
      left_join(grid, by=join_by(lon_idx, lat_idx)) %>%
      filter(!is.na(ADM2_EN))
    # retrieve data from lon, lat, and time index
    dat[[i]] <- dat[[i]] %>%
      mutate(longitude=lon[lon_idx],
             latitude=lat[lat_idx],
             time=dt[time_idx]) 
    
    # aggregate by the US CDC version of epidemiological year and week
    dat[[i]] <- dat[[i]] %>%
      mutate(epi_year=epiyear(time), 
             epi_week=epiweek(time)) %>%
      group_by(ADM2_EN, 
               epi_year, epi_week) %>%
      summarise(mean=mean(Freq),
                min=min(Freq),
                max=max(Freq),
                sd=sd(Freq)) %>%
      ungroup()
    # rename variables for clarity
    dat[[i]] <- dat[[i]] %>% 
      setNames(c(
        colnames(dat[[i]])[1:3],
        paste(vars[i], 
              c("mean", "min", "max", "sd"), sep="_")
      ))
  }
  
  # close the open netCDF file
  nc_close(nc_data)
  
  # combine data frames for different variables using a full join
  dat <- dat %>% 
    reduce(full_join, 
           by = join_by(ADM2_EN, 
                        epi_year, epi_week))
  return(dat)
}

# retrieve land covariate data for the years 2013 to 2022
land_covars <- 
  lapply(2013:2023, land_data_fun)

# garbage collection to reduce memory usage
gc(verbose=TRUE, full=TRUE)

# combine data frames for different years
land_covars <- land_covars %>% 
  reduce(full_join)

# save the extracted land covariates
saveRDS(land_covars, file="spat_temp_covars.rds")
