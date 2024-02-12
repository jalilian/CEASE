
# S2 Prototype Land Cover 20m Map of Africa 2016
# Source: https://2016africalandcover20m.esrin.esa.int/
# Provided to the public by the ESA Climate Change Initiative and 
#      its Land Cover project as the source of the CCI-LC database
# Copyright Notice:
# © Contains modified Copernicus data (2015/2016)
# © ESA Climate Change Initiative - Land Cover project 2017

# format: GeoTIFF

# value codes:
# 0   No data
# 1   Tree cover areas
# 2   Shrubs cover areas
# 3   Grassland
# 4   Cropland
# 5   Vegetation aquatic or regularly flooded
# 6   Lichens Mosses / Sparse vegetation
# 7   Bare areas
# 8   Built-up areas
# 9   Snow and/or Ice
# 10  Open Water

# =========================================================

# Gridded Population of the World (GPW), v4
# resolution: 30 arc second
# provided by Socioeconomic Data and Applications Center, NASA
# https://sedac.ciesin.columbia.edu/data/collection/gpw-v4

# =========================================================

# World Wildlife Fund. Africa: 
# Void-filled digital elevation mode, 2007 [ArcGRID]
# resolution: 15s resolution
# U.S. Geological Survey. 
# Retrieved from https://maps.princeton.edu/catalog/stanford-jm998sr5227 

# =========================================================

get_covars <- local({
  # function to read and extract raster data
  extract_raster_covars <- function(downloaded_file,
                                    data_file,
                                    temp_dir)
  {
    # create a temporary directory to extract the downloaded file
    if (is.null(temp_dir))
      temp_dir <- tempdir()
    if (!dir.exists(temp_dir))
    {
      dir.create(temp_dir)
    }
    
    file_path <- paste0(temp_dir, "/", data_file)
    if (!file.exists(file_path))
    {
      # extract the downloaded file
      unzip(zipfile=downloaded_file, exdir=temp_dir)
    }
    
    # read the raster data using the 'terra' package
    raster <- terra::rast(paste0(temp_dir, "/", data_file))
    
    return(raster)
  }
  
  # function to extract land cover type by coordinates
  files <- list(
    land_cover=c(
      "ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.zip",
      "ESACCI-LC-L4-LC10-Map-20m-P1Y-2016-v1.0.tif",
      "categorical"
    ),
    population_density=c(
      "gpw-v4-population-density-rev11_2020_30_sec_tif.zip",
      "gpw_v4_population_density_rev11_2020_30_sec.tif",
      "numeric"
    ),
    elevation=c(
      "data.zip",
      "af_dem_15s/",
      "numeric"
    )
  )
  
  
  get_land_covers_points <- function(coords, 
                                     path,
                                     temp_dir=NULL)
  {
    covars <- vector("list", length=length(files))
    for (j in 1:length(files))
    {
      r <- extract_raster_covars(paste0(path, files[[j]][1]),
                                 files[[j]][2],
                                 temp_dir)
      covars[[j]] <- terra::extract(r, coords)
    }
    data.frame(covars)
  }
  
  get_land_covars_polygon <- function(map,
                                      path,
                                      temp_dir=NULL)
  {
    covars <- vector("list", length=length(files))
    for (j in 1:length(files))
    {
      r <- extract_raster_covars(paste0(path, files[[j]][1]),
                                 files[[j]][2],
                                 temp_dir)
      #map <- sf::st_transform(map, crs(r))
      out <- vector("list", length=nrow(map))
      for (i in 1:nrow(map))
      {
        cat(i, ": ")
        cp <- terra::crop(r, map[i, ], mask=FALSE, touches=TRUE)
        cp <- terra::mask(cp, map[i, ])
        cp <- terra::values(cp)
        cp <- cp[!is.na(cp)]
        
        out[[i]] <- switch(files[[j]][3], numeric={
          c(mean=mean(cp), min=min(cp), max=max(cp), sd=sd(cp))
        }, categorical={
          prop.table(table(cp))
        })
        cat("done.\n")
      }
      out <- data.frame(out)
      if (files[[j]][3] == "categorical")
        out[is.na(out)] <- 0
    }
    data.frame(covars)
  }
  
  get_land_covars <- function(input, path, temp_dir=NULL)
  {
    if (is(input, "matrix"))
    {
      get_land_covers_points(input, path, temp_dir)
    } else{
      if (is(input, "sf"))
        get_land_covars_polygon(input, path, temp_dir)
      else
        print("wrong input type")
    }
  }
})

# Examples
if (FALSE)
{
  get_covars(cbind(runif(10, -10, 10), runif(10, 10, 20)),
             path="~/Downloads/Africa_covars/")
  
  setwd(tempdir())
  download.file(
    "https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_ETH_shp.zip",
    destfile="map.zip"
    )
  unzip("map.zip", exdir="map/")
  library("sf")
  eth_map <- read_sf("map/", layer="gadm41_ETH_1")
  get_covars(eth_map,
             path="~/Downloads/Africa_covars/")
}
