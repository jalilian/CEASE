
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

# =========================================================
# spatiotemporal covariates
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
  
  grid <- 
    expand_grid(lon_idx=1:length(lon),
                lat_idx=1:length(lat)) %>%
    rowwise() %>% 
    mutate(points=list(st_point(c(lon[lon_idx], 
                                  lat[lat_idx])))) %>%
    st_as_sf()
  
  st_crs(grid) <- st_crs(eth_map)
  
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
    
    dimnames(vals) <- 
      list(lon_idx=1:length(lon), 
           lat_idx=1:length(lat),
           time_idx=1:length(dt))
    
    dat[[i]] <- as.data.frame.table(vals) %>%
      mutate(lon_idx=as.integer(lon_idx),
             lat_idx=as.integer(lat_idx),
             time_idx=as.integer(time_idx))
    
    # determine grid points inside map regions
    dat[[i]] <- dat[[i]] %>%
      left_join(grid, by=join_by(lon_idx, lat_idx)) %>%
      filter(!is.na(ADM2_EN))
    # retrieve data from lon, lat, and time index
    dat[[i]] <- dat[[i]] %>%
      mutate(longitude=lon[lon_idx],
             latitude=lat[lat_idx],
             time=dt[time_idx]) 
    
    # aggregate by map regions and 
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
    # rename columns
    dat[[i]] <- dat[[i]] %>% 
      setNames(c(
        colnames(dat[[i]])[1:3],
        paste(vars[i], 
              c("mean", "min", "max", "sd"), sep="_")
      ))
  }
  
  # close the open netCDF file
  nc_close(nc_data)
  
  dat <- dat %>% 
    reduce(full_join, 
           by = join_by(ADM2_EN, 
                        epi_year, epi_week))
  return(dat)
}

land_covars <- 
  lapply(2013:2022, land_data_fun)

# garbage collection to reduce memory usage
gc(verbose=TRUE, full=TRUE)

land_covars <- land_covars %>% 
  reduce(full_join)
