
# =========================================================
# European Centre for Medium-Range Weather Forecasts (ECMWF) 
# Copernicus's Climate Data Store (CDS)
# https://cds.climate.copernicus.eu
# DOI: 10.24381/cds.e2161bac

# spatial resolution: 0.1x0.1 degrees
# temporal resolution: monthly
# temporal availability: from January 1981 to 2-3 months before the present
# =========================================================

get_csd <- local({
  library("tidyverse")
  library("ecmwfr")
  library("ncdf4")
  library("ncdf4.helpers")
  
  # climate variables 
  cvars <- c("10m_u_component_of_wind", 
             "10m_v_component_of_wind", 
             "leaf_area_index_high_vegetation", 
             "leaf_area_index_low_vegetation", 
             "skin_temperature", 
             "surface_pressure",
             "total_precipitation", 
             "volumetric_soil_water_layer_1")
  
  get_data <- function(user, cds.key,
                       year, month, area, 
                       temp_dir=NULL)
  {
    # set secret ECMWF token
    wf_set_key(user=user, key=cds.key, service="cds")
    
    # create a temporary directory to extract the downloaded file
    if (is.null(temp_dir))
      temp_dir <- tempdir()
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
      variable = cvars,
      # temporal framework: year, month, day, hour
      year = as.character(year),
      month = as.character(month),#sprintf("%02d", 1:12),
      day = sprintf("%02d", 1:31),
      time = "12:00",
      # geographical region
      #      North, West, South, East
      area = area,
      # output file format
      format = "netcdf.zip",
      # output file name
      target = paste0("cds_hourly_", year, ".zip")
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
    unzip(zipfile=paste0(temp_dir, "cds_hourly_", year, ".zip"), 
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
      # retrieve data from lon, lat, and time index
      dat[[i]] <- dat[[i]] %>%
        mutate(longitude=lon[lon_idx],
               latitude=lat[lat_idx],
               time=dt[time_idx]) 
      
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
  
})

if (FALSE)
{
  # user credentials for ECMWF data access
  user <- "****************"
  cds.key <- "********************************"
  
  #      North, West, South, East
  area <- c(10, 30, 5, 35)
  a <- get_csd(user, cds.key, 2020, 8, area)
}
