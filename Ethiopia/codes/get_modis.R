
get_modis <-  local({
  
  library("tidyverse")
  library("rstac")
  library("terra")
  
  get_modis_bbox <- function(collections, asset_key,
                             bbox, datetime, crs=NULL, 
                             output_dir=tempdir())
  {
    # STAC search API
    search_results <- 
      stac("https://planetarycomputer.microsoft.com/api/stac/v1") %>%
      stac_search(
        # collection IDs to include in the search for items
        collections = collections,
        # bounding box (xmin, ymin, xmax, ymax)
        bbox = bbox, 
        # date-time range
        datetime = datetime, 
        # maximum number of results
        limit = 999
      ) %>%
      # HTTP GET requests to STAC web services
      get_request() %>%
      # allow access assets from Microsoft's Planetary Computer
      items_sign(sign_fn=sign_planetary_computer())
    
    # fetch items all STAC Items
    items <- items_fetch(search_results)
    # extract asset IDs
    asset_ids <- sapply(items$features, function(x) { x$id })
    asset_ids <- strsplit(asset_ids, split="\\.")
    
    # dates
    dates <- as.Date(sapply(asset_ids, function(x) x[2]), 
                     format="A%Y%j")
    # extract year and month
    dates <- factor(paste(year(dates),  
                       month(dates, label=FALSE), 
                       sep="-"))
    # horizontal tile number, vertical tile number
    # tiles are 10 degrees by 10 degrees at the equator
    hv <- factor(sapply(asset_ids, function(x) x[3]))

    # download items
    download_items <- items %>%
      assets_download(assets_name="map", 
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
    assets <- tibble(hv, dates, asset_links) %>%
      mutate(across(all_of(asset_key), ~ map(.x, rast)))
    
    assets %>%
      group_by(hv, dates) %>% 
      summarize(across(all_of(asset_key), 
                       ~ list(app(rast(.x), mean, na.rm=TRUE))
                       ),
                .groups="drop") %>%
      group_by(dates) %>% 
      summarize(across(all_of(asset_key),
                       ~ list(Reduce(function(x, y) mosaic(x, y, fun="mean"), 
                                     .x))
                       ))
  }
  
  get_modis_points <- function(collections, asset_key,
                               coords, crs="EPSG:4326", 
                               datetime, output_dir=tempdir())
  {
    coords <- vect(xy, crs=crs)
    bbox <- as.vector(terra::ext(coords))
    bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])

    get_modis_bbox(collections, asset_key, bbox, datetime,
                   output_dir=output_dir) %>% 
    group_by(dates) %>% 
      summarize(across(all_of(asset_key),
                       ~ lapply(.x, function(o){ 
                         terra::extract(terra::project(o,  crs), 
                                        coords, ID=FALSE, xy=TRUE)
                       })
      )) %>%
      unnest(cols=all_of(asset_key), names_sep="__") %>%
      rename_with(~ gsub("__mean", "", .x), 
                  .cols = contains("__mean")) %>%
      select(-matches("__x$") | matches("__x$")[1]) %>%
      rename(long = matches("__x$")[1]) %>%
      select(-matches("__y$") | matches("__y$")[1]) %>%
      rename(lat = matches("__y$")[1]) %>%
      relocate(long, lat, .after=dates)
  }
  
  get_modis_map <- function(collections, asset_key,
                            map, crs="EPSG:4326", 
                            datetime, output_dir=tempdir())
  {
    library("sf")
    cmap <- st_transform(map, crs=crs)
    bbox <- st_bbox(cmap)
    bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])
    
    
    get_modis_bbox(collections, asset_key, bbox, datetime,
                          output_dir=output_dir) %>% group_by(dates) %>% 
      summarize(across(all_of(asset_key),
                       ~ lapply(.x, function(o){ 
                         terra::crop(terra::project(o,  crs), cmap)
                       })
      ))
  }
  
  get_modis_data <- function(collections, asset_key,
                             what, crs="EPSG:4326", 
                             datetime, output_dir=tempdir())
  {
    switch(class(what)[1], numeric={
      if (length(what) == 4)
        get_modis_bbox(collections=collections, 
                       asset_key=asset_key,
                       bbox=what, datetime=datetime, 
                       crs=crs, output_dir=output_dir)
      else
        stop("numeric vector of length 4 is required")
    }, matrix={
      if (ncol(what) != 2)
        stop("a matrix of coordinates with two columns (Long, Lat) is required")
      get_modis_points(collections=collections, 
                       asset_key=asset_key,
                       coords=what, crs=crs, 
                       datetime=datetime, 
                       output_dir=output_dir)
    }, sf={
      get_modis_map(collections=collections, 
                    asset_key=asset_key,
                    map=what, crs=crs, 
                    datetime=datetime, 
                    output_dir=output_dir)
    })
  }
})

if (FALSE)
{
  # MODIS Land Surface Temperature/Emissivity 8-Day
  aa <- get_modis(collections="modis-11A2-061", 
                  asset_key=c("LST_Day_1km", # 8-day daytime 1km grid Landsurface Temperature
                              "LST_Night_1km" # 8-day nighttime 1km grid Landsurface Temperature
                  ), 
                  what=c(30, 2, 32, 4),
                  datetime="2023-11-01/2024-02-28")
  
  # MODIS Gross Primary Productivity 8-Day Gap-Filled
  bb <- get_modis(collections="modis-17A2HGF-061", 
                  asset_key=c("Gpp_500m", #Gross Primary Productivity
                              "PsnNet_500m" # Net Photosynthesis
                  ), 
                  what=c(30, 2, 32, 4),
                  datetime="2023-11-01/2024-02-28")
  plot(bb$Gpp_500m[[1]], col = rev(terrain.colors(100)), colNA="grey")
  
  xy <- cbind(runif(100, 30, 32), runif(100, 2, 4))
  # MODIS Net Evapotranspiration Yearly Gap-Filled
  cc <- get_modis(collections="modis-16A3GF-061", 
                  asset_key=c("ET_500m", # Total of Evapotranspiration
                              "LE_500m", # Average of Latent Heat Flux
                              "PET_500m" # Total Potential Evapotranspiration
                  ),
                  what=xy,
                  datetime="2023-11-01/2024-02-28")
  
  # MODIS Burned Area Monthly
  map <- sf::read_sf("https://geodata.ucdavis.edu/gadm/gadm4.1/kmz/gadm41_IRN_1.kmz")
  map <- map[c(14, 18, 10, 12, 17), ]
  dd <- get_modis(collections="modis-64A1-061",
                  asset_key = c("Burn_Date", # Burn day of year
                                "Burn_Date_Uncertainty"# Estimated uncertainty in burn day
                  ),
                  what=map,
                  datetime="2024-06-01/2024-09-30")
  plot(dd$Burn_Date[[1]], col=rev(heat.colors(100)), 
       colNA = "black", main="June 2024")
  plot(st_geometry(map), add=TRUE, border="blue")
}
