
library("sf")
library("tidyverse")
library("readxl")

# =========================================================
# Ethiopian Public Health Institute (EPHI) 
# weekly malaria surveillance data from 2013 to 2022
# restricted access data

# path to the malaria data file
data_path <- "~/Downloads/Ethiopia/"

# reading the excel file containing the data
eth_data <- 
  read_excel(paste0(data_path,
                    "EPHI_Malaria Surveillance Data from 2013-2022.xlsx"), 
             sheet="Sheet1",
             range=cell_limits(c(1, 1), 
                               c(342472, 13)), 
             col_names = TRUE,
             na = "") %>%
  # transform Month to a factor
  mutate(Month=factor(Month, levels=month.name)) %>%
  # convert year and week to date
  #!!! what does Epidemic week means?
  # Is it week starting from 1st of January each year?
  mutate(date1=ymd(paste(Year, "01", "01", sep="-")) + 
           weeks(Epidemic_Week)) %>%
  # Or is it simply week of the year? (more likely)
  mutate(date2=parse_date_time(paste(Year, Epidemic_Week, 1, sep="-"),
                          "Y-W-w"))
#!!! Issue 1: Year 2015 has 53 weeks instead of 52 weeks
#!!! cause 862 NA in date for year 2015
#!!! needs to be fixed

#!!! Issue 2: some areas (Regions/Zones/Woredas) have less than 52 week data
#!!! Are these missing records? Or they are aggregated with other weeks?


#!!! Issue 3: No data is available from 2020 to 2022
#!!! for 2022, only 4 weeks in January , 3 weeks in July and 3 weeks in August
#!!! Is this due to COVID-19?

# Extract 2013-2019: 
# excluding missing weeks in 2020-2021 and spars weeks in 2022
eth_data <- eth_data %>%
  filter(Year <= 2019) 

# =========================================================
# read map of administrative divisions of Ethiopia in 2015
# data from a shapefile provided by the Stanford Digital Repository

eth_map <- readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2015_Stanford.rds"))

# =========================================================
# Matching Zones in the malaria survilance dataset with the 2015 map

eth_data <- eth_data %>%
  # exclude records with ambiguous region name: Regional
  filter(ZoneName != "Regional") %>%
  # change spelling of some regions
  mutate(RegionName=
           case_match(RegionName,
                      "Gambella" ~ "Gambela",
                      "SNNPR" ~ "SNNP",
                      .default=RegionName)) %>%
  # change name and spelling of some zones
  # stor modified zone names in the new variable ZoneName2
  mutate(ZoneName2=
           # region by region
           case_match(ZoneName,
                      # Addis Ababa region
                      "Addis Ketema" ~ "Addis Ababa",
                      "Akaki Kaliti" ~ "Akaki Kaliti",
                      "Arada" ~ "Addis Ababa",
                      "Bole" ~ "Addis Ababa",
                      "Chirkos" ~ "Addis Ababa",
                      "Gulele" ~ "Addis Ababa",
                      "Kolfe Keraniyo" ~ "Addis Ababa",
                      "Lideta" ~ "Addis Ababa",
                      "Nefas Silk Lafto" ~ "Addis Ababa",
                      "Yeka" ~ "Addis Ababa",
                      # Afar region
                      "Zone 01" ~ "Afar Zone 1",
                      "Zone 02" ~ "Afar Zone 2",
                      "Zone 03" ~ "Afar Zone 3",
                      "Zone 04" ~ "Afar Zone 4",
                      "Zone 05" ~ "Afar Zone 5",
                      # Amhara region
                      "Awi" ~ "Agew Awi",
                      "Argoba Special Woreda" ~ "Argoba",
                      "Bahir Dar Liyu Town" ~ "Bahir Dar Special Zone",
                      "Dese Town" ~ "South Wollo",
                      "Centeral  Gondar" ~ "North Gondar",
                      "West Gondar" ~ "North Gondar",
                      "Gonder Town" ~ "North Gondar",
                      "Oromiya" ~ "Oromia",
                      # Benishangul-Gumuz region
                      "Assosa" ~ "Asosa",
                      "Maokomo Special" ~ "Asosa",
                      # Dire Dawa region
                      "Dredewa" ~ "Dire Dawa",
                      # Gambela region
                      "Agnuwak" ~ "Agnuak",
                      "Mejenger" ~ "Majang",
                      "Nuwer" ~ "Nuer",
                      "Etang Spe." ~ "Agnuak",
                      # Oromia region
                      "Adama Special Town" ~ "East Shewa",
                      "Ambo town" ~ "West Shewa",
                      "Assela Town" ~ "Arsi",
                      "Bishoftu Town" ~ "East Shewa",
                      "Buno Bedele" ~ "Ilubabor",
                      "Burayu Town" ~ "West Shewa",
                      "Dukem Town" ~ "East Shewa",
                      #!!!    "Finfine Zuria" ~ ??,
                      #!!! Finfine is established by merging parts of
                      #!!! several zones around Addis Ababa
                      #!!! fixed in the following by wored names
                      "Gelan Town" ~ "East Shewa",
                      "Holeta town" ~ "West Shewa",
                      "Jimma Spe Town" ~ "Jimma",
                      "Lege Dadi Lege Tafo Town" ~ "North Shewa",
                      "Modjo town" ~ "East Shewa",
                      "Nekemte Town" ~ "East Wellega",
                      "Robe town" ~ "Bale",
                      "Sebeta Town" ~ "South West Shewa",
                      "Shashamane Town" ~ "West Arsi",
                      "Sululta Town" ~ "North Shewa",
                      "Woliso town" ~ "South West Shewa",
                      "Qeleme Wellega" ~ "Kelem Wellega",
                      "west Guji" ~ "Guji",
                      "Ilu Aba Bora" ~ "Ilubabor",
                      # Somali region
                      "FAAFAN" ~ "Fafan",
                      "SITTI" ~ "Sitti",
                      "SHABEELE" ~ "Shabelle",
                      "NOGOB" ~ "Nogob",
                      "Degehabur" ~ "Jarar",
                      "Dhewa" ~ "Liben",
                      "Erar" ~ "Nogob",
                      "Fik" ~ "Nogob",
                      "Gode" ~ "Shabelle", 
                      "Jijiga" ~ "Fafan",
                      "Doollo" ~ "Doolo",
                      "Warder" ~ "Doolo",
                      "Shinile" ~ "Sitti",
                      # SNNP region
                      "Basketo Town" ~ "Basketo",
                      "Dawuro" ~ "Dawro",
                      "Hawassa Town" ~ "Sidama",
                      "Halaba" ~ "Alaba",
                      "Kefa" ~ "Keffa",
                      "Konta Town" ~ "Konta",
                      "Segen" ~ "Konso",
                      "Siliti" ~ "Silti",
                      "Yem Town" ~ "Yem",
                      # Tigray region
                      "Mekele Especial Zone" ~ "South Tigray",
                      "South East" ~ "South Tigray",
                      #
                      .default=ZoneName)) %>%
  #!!! fixing the issue with Finfine zone 
  mutate(ZoneName2= 
           case_match(ZoneName2,
                      "Akaki" ~ "East Shewa",
                      "Berreh" ~ "North Shewa",
                      "Holeta" ~ "West Shewa",
                      "Mullo" ~ "North Shewa",
                      "Sandefa" ~ "North Shewa",
                      "Sebeta Awas" ~ "South West Shewa",
                      "Sululta Town" ~ "North Shewa",
                      "Walmera" ~ "West Shewa",
                      # make records with ambiguous zone name NA
                      "Zone Level" ~ NA,
                      .default=ZoneName2))

# check the matching 
require("stringdist")
# a function for matching based on dissimilarity distance
sdistmatch <- function(a, b)
{
  a <- sort(a)
  b <- sort(b)
  dd <- stringdist::stringdistmatrix(a, b, method="jaccard")
  out <- matrix(NA, nrow=nrow(dd), ncol=ncol(dd))
  for (i in 1:nrow(out))
  {
    out[i, ] <- b[order(dd[i, ])]
  }
  rownames(out) <- a
  return(out)
}

# matching region by region
discrap <- list()
for (region in unique(eth_map$NAME_1))
{
  a <- unique(eth_data %>% 
                filter(RegionName == region) %>% 
                pull(ZoneName))
  b <- unique(eth_map %>% 
                filter(NAME_1 == region) %>% 
                pull(NAME_2))
  discrap[[region]] <- sdistmatch(a, b)
}

# checking discrepancies
# empty (NULL) list means no discrepancy
lapply(discrap, function(o){ 
  a <- colnames(o)
  b <- o[, 1]
  idx <- a != b
  cbind(a[idx], b[idx])
  })
