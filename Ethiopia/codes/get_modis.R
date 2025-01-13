
library("tidyverse")
library("rstac")
library("terra")

get_modis <- local({
  # MODIS collections and variables
  ca <- list(
    # MODIS Land Surface Temperature/Emissivity 8-Day (1km)
    "modis-11A2-061" = c(
      "LST_Day_1km", # 8-day daytime 1km grid Landsurface Temperature
      "LST_Night_1km" # 8-day nighttime 1km grid Landsurface Temperature
    ),
    # MODIS Gross Primary Productivity 8-Day (500m)
    "modis-17A2H-061" = c(
      "Gpp_500m", # Gross Primary Productivity
      "PsnNet_500m" # Net Photosynthesis
    ),
    # MODIS Surface Reflectance 8-Day (500m)
    "modis-09A1-061" = c(
      "sur_refl_b07" # Surface Reflectance Band 7 (2105-2155 nm)
    ),
    # MODIS Net Evapotranspiration Yearly Gap-Filled (500m)
    "modis-16A3GF-061" = c(
      "ET_500m", # Total of Evapotranspiration
      "LE_500m", # Average of Latent Heat Flux
      "PET_500m" # Total Potential Evapotranspiration
    ),
    # MODIS Leaf Area Index/FPAR 8-Day (500m)
    "modis-15A2H-061" = c(
      "Lai_500m", # Leaf Area Index
      "Fpar_500m" # Fraction of Photosynthetically Active Radiation
    ),
    # MODIS Vegetation Indices 16-Day (250m)
    "modis-13Q1-061" = c(
      "250m_16_days_EVI", # 16 day EVI
      "250m_16_days_NDVI" # 16 day NDVI
    ), 
    # MODIS Thermal Anomalies/Fire 8-Day (1km)
    "modis-14A2-061" = c(
      "FireMask" # Confidence of fire
    )
  )
  
  handel_ids <- function(items)
  {
    map2(
      map(items$features, ~ .x$collection),
      map(items$features, ~ .x$id),
      function(collection, id_name) {
        if (startsWith(collection, "modis-")) 
        { 
          # Moderate Resolution Imaging Spectroradiometer (MODIS)
          # Split by "." and extract date and tile
          parts <- strsplit(id_name, "\\.")[[1]]
          # date
          date <- as.Date(parts[2], format = "A%Y%j")
          # horizontal tile number, vertical tile number
          # tiles are 10 degrees by 10 degrees at the equator
          tile <- parts[3]
          return(list(date=date, tile=tile))
        } else if (startsWith(collection, "io-lulc")) 
        {
          # Land Use and Land Cover (LULC) by  Impact Observatory (IO)
          # Split by "-" and extract date and tile
          parts <- strsplit(id_name, "-")[[1]]
          # date
          date <- parts[2]
          # tile
          tile <- parts[1]
          return(list(date=date, tile=tile))
        } else {
          stop("Collection type not recognized")
        }
      }
    ) %>% bind_rows()
    
  }
  
  get_modis_bbox <- function(collections, asset_key,
                             bbox, crs="EPSG:4326", 
                             datetime, crop=TRUE, 
                             aggregate=FALSE,
                             output_dir=tempdir())
  {
    if (!any(startsWith(collections, c("modis-", "io-lulc"))))
      stop("Implemented for the modis and io-lulc products")
      
    # STAC web service: Microsoft Planetary Computer 
    items <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") %>%
      # STAC search API
      stac_search(
        # collection IDs to include in the search for items
        collections = collections,
        # bounding box (xmin, ymin, xmax, ymax) in  WGS84 longitude/latitude
        bbox = bbox, 
        # date-time range
        datetime = datetime, 
        # maximum number of results
        limit = 999
      ) %>%
      # HTTP GET requests to STAC web services
      get_request() %>%
      # allow access assets from Microsoft's Planetary Computer
      items_sign(sign_fn=sign_planetary_computer()) %>%
      # fetch all STAC Items
      items_fetch()
    
    if (!all(unlist(map(items$features, ~ asset_key %in% names(.x$assets)))))
      stop("asset_key is not available for some or all of features")
    
    # select the assets in asset_key
    items <- items %>% 
      assets_select(asset_names=asset_key)

    # extract asset IDs
    asset_ids <- handel_ids(items)
    
    # download items
    download_items <- items %>%
      assets_download(asset_names = asset_key, 
                      items_max=length(items$features), 
                      overwrite=TRUE, output_dir=output_dir)
    
    # extract asset links
    asset_fun <- function(download_items, asset_key)
    {
      sapply(asset_key, function(key){
        sapply(download_items$features, 
               function(feature) {
                 feature$assets[[key]]$href
               })
      }) %>% as_tibble()
    }
    
    asset_links <- asset_fun(download_items, asset_key=asset_key)
    
    # create a tibble with raster data
    assets <- tibble(asset_ids, asset_links) %>%
      mutate(across(all_of(asset_key), ~ map(.x, rast)))
    
    
    # sinusoidal projection CRS used by MODIS
    sincrs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"
    # project the bounding box to the sinusoidal CRS
    sbbox <- rast(terra::ext(bbox, xy=TRUE), crs=crs)
    sbbox <- terra::ext(terra::project(sbbox, sincrs))
    assets <- assets %>%
      group_by(date) %>% 
      summarize(across(all_of(asset_key),
                       ~ {
                         # combine adjacent raster tiles
                         mos <- Reduce(function(x, y) terra::mosaic(x, y, 
                                                                    fun="mean"), 
                                       .x)
                         # crop the mosaic to the bounding box if crop is TRUE
                         if (crop)
                           mos <- crop(mos, sbbox)
                         # project the raster to the target CRS
                         mos <- terra::project(mos, crs)
                         list(mos)
                       }),
                .groups = "drop")
    
    unlink(output_dir, recursive=TRUE)
    if (output_dir == tempdir())
      dir.create(tempdir())
    
    if (aggregate)
    {
      assets <- assets %>%
        mutate(date=paste(year(date), month(date), sep="-")) %>%
        group_by(date) %>% 
        summarize(across(all_of(asset_key), 
                         ~ list(app(rast(.x), mean, na.rm=TRUE))
        ),
        .groups="drop")
    }
    return(assets)
  }
  
  get_modis_points <- function(collections, asset_key,
                               coords, crs="EPSG:4326", 
                               aggregate=FALSE,
                               datetime, output_dir=tempdir())
  {
    # extract the bounding box of coordinates
    bbox <- as.vector(terra::ext(vect(coords, crs=crs)))
    bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])
    # expand the bounding box by 0.1 units in all directions
    bbox <- bbox + c(-1, -1, 1, 1) * 0.1
    
    get_modis_bbox(collections=collections, 
                   asset_key=asset_key, 
                   bbox=bbox, crs=crs, 
                   datetime=datetime,
                   aggregate=aggregate,
                   output_dir=output_dir) %>% 
      group_by(date) %>% 
      summarize(across(all_of(asset_key),
                       ~ lapply(.x, function(o){ 
                         terra::extract(o, coords, 
                                        method="bilinear", 
                                        ID=FALSE, xy=TRUE)
                       })
      )) %>%
      unnest(cols=all_of(asset_key), names_sep="__") %>%
      rename_with(~ gsub("__mean", "", .x), 
                  .cols = contains("__mean")) %>%
      select(-matches("__x$") | matches("__x$")[1]) %>%
      rename(long = matches("__x$")[1]) %>%
      select(-matches("__y$") | matches("__y$")[1]) %>%
      rename(lat = matches("__y$")[1]) %>%
      relocate(long, lat, .after=date)
  }
  
  get_modis_map <- function(collections, asset_key,
                            map, crs="EPSG:4326", 
                            datetime, 
                            aggregate=FALSE,
                            output_dir=tempdir())
  {
    library("sf")
    cmap <- st_transform(map, crs=crs)
    bbox <- st_bbox(cmap)
    bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])
    
    get_modis_bbox(collections=collections, 
                   asset_key=asset_key, 
                   bbox=bbox, 
                   datetime=datetime,
                   aggregate=aggregate,
                   output_dir=output_dir) %>% 
      group_by(date) %>% 
      summarize(across(all_of(asset_key),
                       ~ lapply(.x, function(o){ 
                         terra::crop(o, cmap)
                       })
      ))
  }
  
  
  get_modis_data <- function(collections, asset_key,
                             what, crs="EPSG:4326", 
                             datetime,
                             aggregate=FALSE,
                             output_dir=tempdir())
  {
    switch(class(what)[1], numeric={
      if (length(what) == 4)
        get_modis_bbox(collections=collections, 
                       asset_key=asset_key,
                       bbox=what, datetime=datetime, 
                       crs=crs, aggregate=aggregate,
                       output_dir=output_dir)
      else
        stop("numeric vector of length 4 is required")
    }, matrix={
      if (ncol(what) != 2)
        stop("a matrix of coordinates with two columns (Long, Lat) is required")
      get_modis_points(collections=collections, 
                       asset_key=asset_key,
                       coords=what, crs=crs, 
                       datetime=datetime, 
                       aggregate=aggregate,
                       output_dir=output_dir)
    }, sf={
      get_modis_map(collections=collections, 
                    asset_key=asset_key,
                    map=what, crs=crs, 
                    datetime=datetime, 
                    aggregate=aggregate,
                    output_dir=output_dir)
    })
  }
  
  get_modis_all <- function(collections=names(ca), 
                            asset_key=unname(ca),
                            what, crs="EPSG:4326", 
                            datetime, 
                            aggregate=FALSE,
                            output_dir=tempdir())
  {
    nc <- length(collections)
    if (nc == 1)
    {
      get_modis_data(collections=collections, 
                     asset_key=asset_key,
                     what=what, crs=crs, 
                     datetime=datetime, 
                     aggregate=aggregate,
                     output_dir=output_dir)
    } else{
      if (is.list(asset_key) & (length(asset_key) == nc))
      {
        mapply(get_modis_data, 
               collections=collections,
               asset_key=asset_key,
               MoreArgs=list(what=what, datetime=datetime, aggregate=aggregate),
               SIMPLIFY=FALSE)
      } else{
        stop("collections and asset_keys must have the same length")
      }
    }
  }
})

if (FALSE)
{
  # with  bounding box 
  # MODIS Vegetation Indices 16-Day (250m)
  a1 <- get_modis(collections="modis-13Q1-061", 
                  asset_key=c("250m_16_days_EVI", # 16 day EVI
                              "250m_16_days_NDVI" # 16 day NDVI
                  ), 
                  what=c(-3, 5, -2, 6),
                  datetime="2024-01-01/2024-03-01")
  par(mfrow=c(1, 2))
  plot(a1$`250m_16_days_EVI`[[1]], main=as.character(a1$date[[1]]))
  plot(a1$`250m_16_days_NDVI`[[1]], main=as.character(a1$date[[1]]))
  
  # all selected MODIS variables
  a2 <- get_modis(what=c(-3, 5, -2, 6),
                  datetime="2023-11-01/2024-02-28")
  
  plot(a2$`modis-14A2-061`$FireMask[[1]])
  plot(a2$`modis-09A1-061`$sur_refl_b07[[1]])

  # with coordinates
  xy <- cbind(runif(100, -3, 0), runif(100, 5, 10))
  # MODIS Net Evapotranspiration Yearly Gap-Filled
  a3 <- get_modis(collections="modis-16A3GF-061", 
                  asset_key=c("ET_500m", # Total of Evapotranspiration
                              "LE_500m", # Average of Latent Heat Flux
                              "PET_500m" # Total Potential Evapotranspiration
                  ),
                  what=xy,
                  datetime="2020-01-01/2025-01-01")
  a3
  a3 %>% count(date)

  # with map
  map <- sf::read_sf("https://geodata.ucdavis.edu/gadm/gadm4.1/kmz/gadm41_GHA_0.kmz")
  # MODIS Land Surface Temperature/Emissivity 8-Day
  a4 <- get_modis(collections="modis-11A2-061",
                  asset_key=c("LST_Day_1km", # 8-day daytime 1km grid Landsurface Temperature
                              "LST_Night_1km" # 8-day nighttime 1km grid Landsurface Temperature
                              ),
                  what=map,
                  datetime="2023-12-01/2024-02-01")
  plot(a4$LST_Day_1km[[1]], main=as.character(a4$date[[1]]))
  plot(a4$LST_Night_1km[[1]], main=as.character(a4$date[[1]]))
  
  # Land cover 
  a5 <- get_modis(collections="io-lulc-annual-v02", 
                  asset_key="data", 
                  what=c(-3, 5, -2, 6),
                  datetime=NULL)

  
  # MODIS Burned Area Monthly
  map <- sf::read_sf("https://geodata.ucdavis.edu/gadm/gadm4.1/kmz/gadm41_IRN_1.kmz")
  a6 <- get_modis(collections="modis-64A1-061",
                  asset_key = c("Burn_Date", # Burn day of year
                                "Burn_Date_Uncertainty"# Estimated uncertainty in burn day
                  ),
                  what=map[14, ],
                  datetime="2024-06-01/2024-10-30")
  par(mfrow=c(2, 2))
  plot(a6$Burn_Date[[1]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[1]]))
  plot(st_geometry(map), add=TRUE, border="blue")
  plot(a6$Burn_Date[[2]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[2]]))
  plot(st_geometry(map), add=TRUE, border="blue")
  plot(a6$Burn_Date[[3]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[3]]))
  plot(st_geometry(map), add=TRUE, border="blue")
  plot(a6$Burn_Date[[4]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[4]]))
  plot(st_geometry(map), add=TRUE, border="blue")
  
  a7 <- get_modis(collections="io-lulc-annual-v02", 
                  asset_key="data", 
                  what=c(46, 34, 47, 35),
                  datetime="2023-01-01")
}
