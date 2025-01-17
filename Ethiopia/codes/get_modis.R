
library("tidyverse")
library("rstac")
library("terra")

local({
  
  fill_in_gap <- function(r, w=NULL)
  {
    if (is.null(w))
      w <- 1 / (1 + outer(seq(-4, 4), seq(-4, 4), 
                          function(x, y){ (x^2 + y^2) }))
    while (any(is.na(values(r)))) 
    {
      # apply the focal mean to fill NA values
      r <- focal(r, w=w, fun="mean", na.rm=TRUE, na.policy="only")
    }
    return(r)
  }
  
  handel_ids <- function(items)
  {
    map_dfr(items$features, function(x) 
    {
      collection <- x$collection
      id_name <- x$id
      
      if (startsWith(collection, "modis-"))
      {
        # Moderate Resolution Imaging Spectroradiometer (MODIS)
        # Split by "." and extract date and tile
        parts <- strsplit(id_name, "\\.")[[1]]
        # date: year and day of the year
        date <- as.Date(parts[2], format = "A%Y%j")
        # horizontal tile number, vertical tile number
        # tiles are 10 degrees by 10 degrees at the equator
        tile <- parts[3]
        return(data.frame(date=date, tile=tile))
      } else if (startsWith(collection, "io-lulc")) 
      {
        # Land Use and Land Cover (LULC) by  Impact Observatory (IO)
        # Split by "-" and extract date and tile
        parts <- strsplit(id_name, "-")[[1]]
        # date: year
        date <- parts[2]
        # UTM tile
        # tiles are 6 degrees of longitude wide for north-south zones
        tile <- parts[1]
        return(data.frame(date=date, tile=tile))
      } else {
        stop("Collection type not recognized")
      }
    })
  }
  
  get_plancomp_box <- function(collections, asset_key,
                               bbox, datetime, 
                               crs="EPSG:4326", 
                               crop=TRUE, 
                               fill_in=FALSE,
                               aggregate=FALSE,
                               output_dir=tempdir(),
                               clean_dir=FALSE)
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
                         {
                           sbbox <- rast(terra::ext(bbox, xy=TRUE), crs=crs)
                           sbbox <- terra::project(sbbox, crs(mos, proj=TRUE))
                           mos <- crop(mos, terra::ext(sbbox))
                         }
                         # project the raster to the target CRS
                         mos <- terra::project(mos, crs)
                         list(mos)
                       }),
                .groups = "drop")
    
    if (aggregate)
    {
      assets <- assets %>%
        mutate(date=paste(year(date), month(date), sep="-")) %>%
        group_by(date) %>% 
        summarize(across(all_of(asset_key), 
                         ~ list(app(rast(.x), mean, na.rm=TRUE))
        ), .groups="drop")
    }
    
    if (fill_in)
    {
      assets <- assets %>%
        group_by(date) %>% 
        summarize(across(all_of(asset_key), 
                         ~ lapply(.x, fill_in_gap)
        ), .groups="drop") %>%
        rename_with(~ gsub("__focal", "", .))
    }
    
    if (clean_dir)
    {
      unlink(output_dir, recursive=TRUE)
      if (output_dir == tempdir())
        dir.create(tempdir())      
    }
    
    return(assets)
  }
  
  get_plancomp_points <- function(collections, asset_key,
                                  coords, datetime,
                                  crs="EPSG:4326", 
                                  fill_in=FALSE,
                                  aggregate=FALSE,
                                  output_dir=tempdir(),
                                  clean_dir=FALSE)
  {
    coords <- vect(coords, crs=crs)
    # extract the bounding box of coordinates
    bbox <- as.vector(terra::ext(coords))
    bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])
    # expand the bounding box by 0.1 units in all directions
    bbox <- bbox + c(-1, -1, 1, 1) * 0.2
    
    get_plancomp_box(collections=collections, 
                     asset_key=asset_key, 
                     bbox=bbox, 
                     datetime=datetime,
                     crs=crs,
                     crop=TRUE,
                     fill_in=fill_in,
                     aggregate=aggregate,
                     output_dir=output_dir,
                     clean_dir=clean_dir) %>% 
      group_by(date) %>% 
      summarize(across(all_of(asset_key),
                       ~ lapply(.x, function(o){ 
                         
                         terra::extract(o, coords, 
                                        method="simple", 
                                        ID=FALSE)
                       })
      )) %>%
      mutate(long=lapply(date, function(a) geom(coords)[, "x"]),
             lat=lapply(date, function(a) geom(coords)[, "y"])) %>%
      unnest(cols=all_of(c(asset_key, "long", "lat")),
             names_sep="__") %>%
      rename_with(~ gsub("_mean$", "", .)) %>%
      relocate(long, lat, .after=date)
  }
  
  get_plancomp_map <- function(collections, asset_key,
                               map,  
                               datetime, 
                               crs="EPSG:4326",
                               fill_in=FALSE,
                               aggregate=FALSE,
                               output_dir=tempdir(),
                               clean_dir=FALSE)
  {
    library("sf")
    cmap <- st_transform(map, crs=crs)
    bbox <- st_bbox(cmap)
    bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])
    
    get_plancomp_box(collections=collections, 
                     asset_key=asset_key, 
                     bbox=bbox, 
                     datetime=datetime,
                     crs=crs,
                     crop=TRUE,
                     fill_in=fill_in,
                     aggregate=aggregate,
                     output_dir=output_dir,
                     clean_dir=clean_dir) %>% 
      group_by(date) %>% 
      summarize(across(all_of(asset_key),
                       ~ lapply(.x, function(o){ 
                         terra::crop(o, cmap)
                       })
      ))
  }
  
  
  get_plancomp_data <- function(collections, asset_key,
                                what, 
                                datetime,
                                crs="EPSG:4326",
                                crop=TRUE,
                                fill_in=FALSE,
                                aggregate=FALSE,
                                output_dir=tempdir(),
                                clean_dir=FALSE)
  {
    switch(class(what)[1], numeric={
      if (length(what) == 4)
        get_plancomp_box(collections=collections, 
                         asset_key=asset_key, 
                         bbox=what, 
                         datetime=datetime,
                         crs=crs,
                         crop=crop,
                         fill_in=fill_in,
                         aggregate=aggregate,
                         output_dir=output_dir,
                         clean_dir=clean_dir)
      else
        stop("numeric vector of length 4 is required")
    }, matrix={
      if (ncol(what) != 2)
        stop("a matrix of coordinates with two columns (Long, Lat) is required")
      get_plancomp_points(collections=collections, 
                          asset_key=asset_key, 
                          coords=what, 
                          datetime=datetime,
                          crs=crs,
                          fill_in=fill_in,
                          aggregate=aggregate,
                          output_dir=output_dir,
                          clean_dir=clean_dir)
    }, sf={
      get_plancomp_map(collections=collections, 
                       asset_key=asset_key, 
                       map=what, 
                       datetime=datetime,
                       crs=crs,
                       fill_in=fill_in,
                       aggregate=aggregate,
                       output_dir=output_dir,
                       clean_dir=clean_dir)
    })
  }
  
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
  
  get_plancomp_all <- function(collections=names(ca), 
                               asset_key=unname(ca),
                               what,  
                               datetime,
                               crs="EPSG:4326",
                               crop=TRUE,
                               fill_in=FALSE,
                               aggregate=FALSE,
                               output_dir=tempdir(),
                               clean_dir=FALSE)
  {
    nc <- length(collections)
    if (nc == 1)
    {
      get_plancomp_data(collections=collections, 
                        asset_key=asset_key, 
                        what=what, 
                        datetime=datetime,
                        crs=crs,
                        crop=crop,
                        fill_in=fill_in,
                        aggregate=aggregate,
                        output_dir=output_dir,
                        clean_dir=clean_dir)
    } else{
      if (is.list(asset_key) & (length(asset_key) == nc))
      {
        mapply(get_plancomp_data, 
               collections=collections,
               asset_key=asset_key,
               MoreArgs=list(what=what, 
                             datetime=datetime, 
                             crs=crs,
                             crop=crop,
                             fill_in=fill_in, 
                             aggregate=aggregate,
                             output_dir=output_dir,
                             clean_dir=clean_dir),
               SIMPLIFY=FALSE)
      } else{
        stop("collections and asset_keys must have the same length")
      }
    }
  }
  
  get_modis <- get_plancomp_all
  
  # 10m Annual Land Use Land Cover (9-class) V2
  #  by Impact Observatory, Microsoft, and Esri, displays a global map of land use and land cover 
  # classification
  # 1  Water: 
  #     Areas where water was predominantly present throughout the year
  # 2  Trees: 
  #    Areas with significant clustering of tall, dense vegetation
  # 3  Flooded Vegetation: 
  #     Areas of any type of vegetation with obvious intermixing of water throughout a majority of the year
  # 4  Crops: 
  #     Human-planted cereals, grasses, and crops not at tree height
  # 5  Built Area: 
  #     Human-made structures, major road and rail networks, and large homogenous impervious surfaces.
  # 6  Bare Ground: 
  #     Areas of rock or soil with very sparse to no vegetation
  # 7  Snow/Ice: 
  #     Large homogenous areas of permanent snow or ice
  # 8  Clouds: 
  #     Areas with persistent cloud cover, where no land cover information is available.
  # 9  Rangeland: 
  #     Open areas covered in homogenous grasses with little to no taller vegetation.

  
  get_lulc <- function(what, datetime, 
                       crs="EPSG:4326",
                       fact=0,
                       output_dir=tempdir(),
                       clean_dir=FALSE)
  {
    switch(class(what)[1], numeric={
      if (length(what) == 4)
      {
        bbox <- what
      } else
        stop("numeric vector of length 4 is required")
    }, matrix={
      if (ncol(what) != 2)
        stop("a matrix of coordinates with two columns (Long, Lat) is required")
      coords <- vect(what, crs=crs)
      # extract the bounding box of coordinates
      bbox <- as.vector(terra::ext(coords))
      bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])
      # expand the bounding box by 0.1 units in all directions
      bbox <- bbox + c(-1, -1, 1, 1) * 0.25
    }, sf={
      library("sf")
      cmap <- st_transform(what, crs=crs)
      bbox <- st_bbox(cmap)
      bbox <- unname(bbox[c("xmin", "ymin", "xmax", "ymax")])
    },
    stop("Unsupported input type for 'what'")
    )
    
    dat <- get_plancomp_box(collections="io-lulc-annual-v02", 
                            asset_key="data",
                            bbox=bbox, 
                            datetime=datetime,
                            crs=crs,
                            crop=TRUE,
                            fill_in=FALSE,
                            aggregate=FALSE,
                            output_dir=output_dir,
                            clean_dir=clean_dir)
    
    # 0   No Data
    # 1   Water
    # 2   Trees
    # 4   Flooded vegetation
    # 5   Crops
    # 7   Built area
    # 8   Bare ground
    # 9   Snow/ice
    # 10  Clouds
    # 11  Rangeland
    # land cover classification levels
    levs <- c("Water", "Trees", NA, "Flooded vegetation", 
              "Crops", NA, "Built area", "Bare ground", 
              "Snow/ice", "Clouds", "Rangeland")
    
    dat <- dat %>%
      mutate(data = lapply(data, function(r) 
      {
        v <- values(r)
        # non-finite or integer values to NA
        v[(!is.finite(v)) | (v %% 1 != 0)] <- NA
        values(r) <- levs[v]
        if (fact > 0)
          r <- aggregate(r, fact=fact, fun="modal", na.rm=TRUE)
        return(r)
      }))
    
    if (class(what)[1] == "matrix")
    {
      dat <- dat %>% 
        group_by(date) %>% 
        summarize(across(all_of("data"),
                         ~ lapply(.x, function(o){ 
                           
                           terra::extract(o, coords, 
                                          method="simple", 
                                          ID=FALSE)
                         })
        )) %>%
        mutate(long=lapply(date, function(a) geom(coords)[, "x"]),
               lat=lapply(date, function(a) geom(coords)[, "y"])) %>%
        unnest(cols=c("data", "long", "lat")) %>%
        relocate(long, lat, .after=date)
    } else if (class(what)[1] == "sf")
    {
      dat <- dat %>% 
        group_by(date) %>% 
        summarize(across(all_of("data"),
                         ~ lapply(.x, function(o){ 
                           terra::crop(o, cmap)
                         })
        ))
    }
    return(dat)
  }
  
  assign("get_modis", get_modis, envir = .GlobalEnv)
  assign("get_lulc", get_lulc, envir = .GlobalEnv)
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
  a5 <- get_lulc(what=c(-3, 5, -2, 6),
                 fact=10,
                 datetime="2023-01-01")
  xy <- cbind(runif(100, -3, -2), runif(100, 5, 6))
  a6 <- get_lulc(what=xy, fact=10,
                 datetime="2023-01-01")
  a7 <- get_lulc(what=map, fact=10,
                 datetime="2023-01-01")
  
  # MODIS Burned Area Monthly
  map <- sf::read_sf("https://geodata.ucdavis.edu/gadm/gadm4.1/kmz/gadm41_IRN_1.kmz")
  a8 <- get_modis(collections="modis-64A1-061",
                  asset_key = c("Burn_Date", # Burn day of year
                                "Burn_Date_Uncertainty"# Estimated uncertainty in burn day
                  ),
                  what=map[14, ],
                  datetime="2024-06-01/2024-10-30")
  par(mfrow=c(2, 2))
  plot(a8$Burn_Date[[1]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[1]]))
  plot(st_geometry(map), add=TRUE, border="blue")
  plot(a8$Burn_Date[[2]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[2]]))
  plot(st_geometry(map), add=TRUE, border="blue")
  plot(a8$Burn_Date[[3]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[3]]))
  plot(st_geometry(map), add=TRUE, border="blue")
  plot(a8$Burn_Date[[4]], col=rev(heat.colors(100)), 
       colNA = "black", main=as.character(a6$date[[4]]))
  plot(st_geometry(map), add=TRUE, border="blue")

  a9 <- get_lulc(what=c(46.9, 34.2, 47.2, 34.5),
                  datetime="2023-01-01")
  
}
