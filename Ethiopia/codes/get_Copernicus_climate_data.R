
# =========================================================
# European Centre for Medium-Range Weather Forecasts (ECMWF) 
# Copernicus's Climate Data Store (CDS)
# https://cds.climate.copernicus.eu
# DOI: 10.24381/cds.e2161bac

# spatial resolution: 0.1x0.1 degrees
# temporal resolution: monthly
# temporal availability: from January 1981 to 2-3 months before the present
# =========================================================

get_cds <- local({
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
  
  get_cds_area <- function(user, cds.key,
                           year, 
                           month=sprintf("%02d", 1:12),
                           day=sprintf("%02d", 1:31), 
                           time="12:00",
                           area, 
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
      month = as.character(month),
      day = as.character(day),
      time = as.character(time),
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
    unzip(zipfile=paste0("cds_hourly_", year, ".zip"), 
          exdir=paste0(year, "/"), 
          overwrite=TRUE)
    
    # open netCDF file containing the data
    nc_data <- nc_open(paste0(year, "/data.nc"))
    
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
      
      if(length(dim(vals)) < 3)
      {
        vals <- array(vals, 
                      dim=c(length(lon), 
                            length(lat), 
                            length(dt)))
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
               time=dt[time_idx]) %>%
        select(-c(lon_idx, lat_idx, time_idx)) %>%
        relocate(Freq, .after=time) %>%
        # rename variable for clarity
        rename(!!quo_name(vars[i]) := "Freq")
    }
    
    # close the open netCDF file
    nc_close(nc_data)
    
    # combine data frames for different variables using a full join
    dat <- dat %>% 
      reduce(full_join, 
             by = join_by(longitude, latitude, time))
    return(dat)
  }
  
  get_cds_points <- function(user, cds.key,
                             year,
                             month=sprintf("%02d", 1:12),
                             day=sprintf("%02d", 1:31),
                             time="12:00",
                             coords,
                             temp_dir=NULL)
  {
    # coordinates of desired locations
    x1 <- coords[, 1] # longitude
    y1 <- coords[, 2] # latitude
    
    #         North, West, South, East
    area <- c(max(y1), min(x1), min(y1), max(x1))
    
    cds_dat <- get_cds_area(user=user, cds.key=cds.key, 
                            year=year, month=month, 
                            day=day, time=time,
                            area=area, temp_dir=temp_dir)
    
    xy2 <- cds_dat %>% distinct(longitude, latitude)
    x2 <- xy2 %>% pull(longitude)
    y2 <- xy2 %>% pull(latitude)
    # squared Euclidean distance for all desired locations to cds data coordinates
    dd <- outer(x1, x2, "-")^2 + outer(y1, y2, "-")^2
    
    # index of the closest cds data location for each desired location
    idx <- apply(dd, 1, which.min)
    out <- vector("list", length=nrow(coords))
    for (i in 1:length(out))
    {
      out[[i]] <- cds_dat %>%
        filter(longitude == x2[idx[i]] & latitude == y2[idx[i]]) %>%
        mutate(x1=x1, y1=y1) %>%
        relocate(x1, .before=longitude) %>%
        relocate(y1, .before=longitude)
    }

    return(bind_rows(out))
  }
  
  get_cds_map <- function(user, cds.key,
                          year,
                          month=sprintf("%02d", 1:12),
                          day=sprintf("%02d", 1:31),
                          time="12:00",
                          map,
                          temp_dir=NULL)
  {
    library("sf")
    area <- st_bbox(map)
    #         North, West, South, East
    area <- c(area[4], area[1:3])
    
    cds_dat <- get_cds_area(user=user, cds.key=cds.key, 
                            year=year, month=month, 
                            day=day, time=time,
                            area=area, temp_dir=temp_dir)
    
    cds_dat <- st_as_sf(cds_dat, 
                        coords=c("longitude", "latitude"),
                        crs=st_crs("wgs84"))
    
    # conduct a spatial join to determine points inside specific map regions
    cds_dat <- st_join(cds_dat, map, join=st_within) %>%
      na.omit()
    
    return(cds_dat)
  }
  
  get_cds_data <- function(user, cds.key,
                           year, 
                           month=sprintf("%02d", 1:12),
                           day=sprintf("%02d", 1:31), 
                           time="12:00",
                           what,
                           temp_dir=NULL)
  {
    switch(class(what)[1], numeric={
      if (length(what) == 4)
        get_cds_area(user=user, cds.key=cds.key, 
                     year=year, month=month, 
                     day=day, time=time,
                     area=what, temp_dir=temp_dir)
      else
        stop("numeric vector of length 4 is required")
    }, matrix={
      get_cds_points(user=user, cds.key=cds.key, 
                     year=year, month=month, 
                     day=day, time=time,
                     coords=what, temp_dir=temp_dir)
    }, sf={
      get_cds_map(user=user, cds.key=cds.key, 
                  year=year, month=month, 
                  day=day, time=time,
                  map=what, temp_dir=temp_dir)
    })
  }
  
})

if (FALSE)
{
  # user credentials for ECMWF data access
  user <- "****************"
  cds.key <- "********************************"
  
  #      North, West, South, East
  area <- c(10, 30, 5, 35)
  a1 <- get_cds(user, cds.key, year=2020, month=8, what=area)
  
  coords <- cbind(seq(30, 35, l=5), 
        seq(5, 10, l=5))
  a2 <- get_cds(user, cds.key, year=2020, month=8, what=coords)
  
  area <- c(9.2, -79.9, 9.1, -79.8)
  a3 <- get_cds(user, cds.key, year=2020, month=8, what=area)
}
