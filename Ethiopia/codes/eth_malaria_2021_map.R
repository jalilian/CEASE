
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
  #!!! what does Epidemic week mean?
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
# read map of administrative divisions of Ethiopia in 2021
# data from a shapefile provided by OCHA
eth_map <- readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))

# =========================================================
# Matching Zones in the malaria survilance dataset with the 2021 map

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
                      # Somali
                      ADM2_EN == "Afder" ~ "Afder",
                      ADM2_EN == "Degehabur" ~ "Jarar",
                      ADM2_EN == "Fik" ~ "Nogob",
                      ADM2_EN == "Gode" ~ "Shabelle",
                      ADM2_EN == "Jijiga" ~ "Fafan",
                      ADM2_EN == "Korahe" ~ "Korahe",
                      ADM2_EN == "Liben" ~ "Liban",
                      ADM2_EN == "Shinile" ~ "Siti",
                      ADM2_EN == "Warder" ~ "Doolo",
                      ADM2_EN == "Doollo" ~ "Doolo",
                      ADM2_EN == "FAAFAN" ~ "Fafan",
                      ADM2_EN == "Jarar" ~ "Jarar",
                      ADM2_EN == "NOGOB" ~ "Nogob",
                      ADM2_EN == "SHABEELE" ~ "Shabelle",
                      ADM2_EN == "SITTI" ~ "Siti",
                      ADM2_EN == "Dhewa" ~ "Daawa",
                      ADM2_EN == "Erar" ~ "Erar",
                      ADM2_EN == "Nogob" ~ "Nogob",
                      # Tigray 
                      ADM2_EN == "Central Tigray" ~ "Central",
                      ADM2_EN == "Eastern Tigray" ~ "Eastern",
                      ADM2_EN == "Mekele Especial Zone" ~ "Mekelle",
                      ADM2_EN == "North Western Tigray" ~ "North Western",
                      ADM2_EN == "South East" ~ "South Eastern",
                      ADM2_EN == "South Tigray" ~ "Southern",
                      ADM2_EN == "Western Tigray" ~ "Western",
                      # Oromia 
                      ADM2_EN == "Adama Special Town" ~ "East Shewa",
                      ADM2_EN == "Arsi" ~ "Arsi",
                      ADM2_EN == "Assela Town" ~ "Arsi",
                      ADM2_EN == "Bale" ~ "Bale",
                      ADM2_EN == "Bishoftu Town" ~ "East Shewa",
                      ADM2_EN == "Borena" ~ "Borena",
                      ADM2_EN == "Burayu Town" ~ "West Shewa",
                      ADM2_EN == "Dukem Town" ~ "East Shewa",
                      ADM2_EN == "East Hararge" ~ "East Hararge",
                      ADM2_EN == "East Shewa" ~ "East Shewa",
                      ADM2_EN == "East Wellega" ~ "East Wellega",
                      ADM2_EN == "Finfine Zuria" ~ "Finfine Special",
                      ADM2_EN == "Gelan Town" ~ "Finfine Special",
                      ADM2_EN == "Guji" ~ "Guji",
                      ADM2_EN == "Horo Gudru Wellega" ~ "Horo Gudru Wellega",
                      ADM2_EN == "Ilu Aba Bora" ~ "Ilu Aba Bora",
                      ADM2_EN == "Jimma" ~ "Jimma",
                      ADM2_EN == "Jimma Spe Town" ~ "Jimma",
                      ADM2_EN == "Lege Dadi Lege Tafo Town" ~ "Finfine Special",
                      ADM2_EN == "Nekemte Town" ~ "East Wellega",
                      ADM2_EN == "North Shewa" ~ "North Shewa (OR)",
                      ADM2_EN == "Qeleme Wellega" ~ "Kelem Wellega",
                      ADM2_EN == "Sebeta Town" ~ "Finfine Special",
                      ADM2_EN == "Shashamane Town" ~ "West Arsi",
                      ADM2_EN == "South West Shewa" ~ "South West Shewa",
                      ADM2_EN == "Sululta Town" ~ "Finfine Special",
                      ADM2_EN == "West Arsi" ~ "West Arsi",
                      ADM2_EN == "West Hararge" ~ "West Hararge",
                      ADM2_EN == "West Shewa" ~ "West Shewa",
                      ADM2_EN == "West Wellega" ~ "West Wellega",
                      ADM2_EN == "Buno Bedele" ~ "Buno Bedele",
                      ADM2_EN == "west Guji" ~ "West Guj",
                      ADM2_EN == "Ambo town" ~ "West Shewa",
                      ADM2_EN == "Holeta town" ~ "Finfine Special",
                      ADM2_EN == "Modjo town" ~ "East Shewa",
                      ADM2_EN == "Robe town" ~ "Bale",
                      ADM2_EN == "Woliso town" ~ "South West Shewa",
                      # Amhara
                      ADM2_EN == "Argoba Special Woreda" ~ "South Wello",
                      ADM2_EN == "Awi" ~ "Awi",
                      ADM2_EN == "Bahir Dar Liyu Town" ~ "West Gojam",
                      ADM2_EN == "Dese Town" ~ "South Wello",
                      ADM2_EN == "East Gojjam" ~ "East Gojam",
                      ADM2_EN == "Gonder Town" ~ "Central Gondar",
                      ADM2_EN == "North Gondar" ~ "North Gondar",
                      ADM2_EN == "North Shewa" ~ "North Shewa (AM)",
                      ADM2_EN == "North Wollo" ~ "North Wollo",
                      ADM2_EN == "Oromiya" ~ "Oromia",
                      ADM2_EN == "South Gonder" ~ "South Gondar",
                      ADM2_EN == "South Wollo" ~ "South Wello",
                      ADM2_EN == "Wag Himra" ~ "Wag Hamra",
                      ADM2_EN == "West Gojjam" ~ "West Gojam",
                      ADM2_EN == "Centeral  Gondar" ~ "Central Gondar",
                      ADM2_EN == "West Gondar" ~ "West Gondar",
                      # SNNPR
                      ADM2_EN == "Basketo Town" ~ "Basketo",
                      ADM2_EN == "Bench Maji" ~ "Bench Sheko",
                      ADM2_EN == "Dawuro" ~ "Dawuro",
                      ADM2_EN == "Gamo Gofa" ~ "Gamo",
                      ADM2_EN == "Gedeo" ~ "Gedeo",
                      ADM2_EN == "Gurage" ~ "Guraghe",
                      ADM2_EN == "Hadiya" ~ "Hadiya",
                      ADM2_EN == "Halaba" ~ "Halaba",
                      ADM2_EN == "Hawassa Town" ~ "Sidama",
                      ADM2_EN == "Kefa" ~ "Kefa",
                      ADM2_EN == "Kembata Tembaro" ~ "Kembata Tembaro",
                      ADM2_EN == "Konta Town" ~ "Konta Special",
                      ADM2_EN == "Segen" ~ "Burji",
                      ADM2_EN == "Sheka" ~ "Sheka",
                      ADM2_EN == "Sidama" ~ "Sidama",
                      ADM2_EN == "Siliti" ~ "Siltie",
                      ADM2_EN == "South Omo" ~ "South Omo",
                      ADM2_EN == "Wolayita" ~ "Wolayita",
                      ADM2_EN == "Yem Town" ~ "Yem Special",
                      ADM2_EN == "Alle" ~ "Alle",
                      ADM2_EN == "Amaro" ~ "Amaro",
                      ADM2_EN == "Basketo" ~ "Basketo",
                      ADM2_EN == "Burji" ~ "Burji",
                      ADM2_EN == "Derashe" ~ "Derashe",
                      ADM2_EN == "Gamo" ~ "Gamo",
                      ADM2_EN == "Konso" ~ "Konso",
                      ADM2_EN == "Silte" ~ "Siltie",
                      ADM2_EN == "Yem" ~ "Yem Special",
                      #SWEPRS 
                      ADM2_EN == "Bench-Sheko" ~ "Bench Sheko",
                      ADM2_EN == "Dawro  ZHD" ~ "Dawuro",
                      ADM2_EN == "Kaffa" ~ "Kefa",
                      ADM2_EN == "KONTA" ~ "Konta Special",
                      ADM2_EN == "SHEKA" ~ "Sheka",
                      ADM2_EN == "West OMO ZHD" ~ "Mirab Omo",
                      # Sidama 
                      ADM2_EN == "Hawassa Town" ~ "Sidama",
                      ADM2_EN == "Hawella Tulla" ~ "Sidama",
                      ADM2_EN == "Sidama  ZHD" ~ "Sidama",
                      # Gambela 
                      ADM2_EN == "Agnuwak" ~ "Agnewak",
                      ADM2_EN == "Etang Spe." ~ "Itang Special woreda",
                      ADM2_EN == "Mejenger" ~ "Majang",
                      ADM2_EN == "Nuwer" ~ "Nuwer",
                      # Harari 
                      ADM2_EN == "Harari" ~ "Harari",
                      # Dire Dawa 
                      ADM2_EN == "Dredewa" ~ "Dire Dawa urban",
                      # Beneshangul gumuz
                      ADM2_EN == "Assosa" ~ "Assosa",
                      ADM2_EN == "Kamashi" ~ "Kamashi",
                      ADM2_EN == "Maokomo Special" ~ "Mao Komo Special",
                      ADM2_EN == "Metekel" ~ "Metekel",
                      ADM2_EN == "SHEKA" ~ "Sheka",
                      # Afar
                      ADM2_EN == "Zone 01" ~ "Awsi /Zone 1",
                      ADM2_EN == "Zone 02" ~ "Kilbati /Zone2",
                      ADM2_EN == "Zone 03" ~ "Gabi /Zone 3",
                      ADM2_EN == "Zone 04" ~ "Fanti /Zone 4",
                      ADM2_EN == "Zone 05" ~ "Hari /Zone 5",
                      # Addis Ababa
                      ADM2_EN == "Addis Ketema" ~ "Region 14",
                      ADM2_EN == "Akaki Kaliti" ~ "Region 14",
                      ADM2_EN == "Arada" ~ "Region 14",
                      ADM2_EN == "Bole" ~ "Region 14",
                      ADM2_EN == "Chirkos" ~ "Region 14",
                      ADM2_EN == "Gulele" ~ "Region 14",
                      ADM2_EN == "Kolfe Keraniyo" ~ "Region 14",
                      ADM2_EN == "Lideta" ~ "Region 14",
                      ADM2_EN == "Nefas Silk Lafto" ~ "Region 14",
                      ADM2_EN == "Yeka" ~ "Region 14",
                      ADM2_EN == "Lemi Kura\r\r\n" ~ "Region 14",
                      ADM2_EN == "Lemi Kura" ~ "Region 14",
                      .default=ZoneName)) %>%
  #!!! fixing the issue with Finfine zone 
  mutate(ZoneName2= 
           case_match(WoredaName,
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

# old and new names of zones
eth_data %>% count(ZoneName, ZoneName2)  %>% print(n=500)

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
                pull(ZoneName2))
  b <- unique(eth_map %>% 
                filter(NAME_1 == region) %>% 
                pull(NAME_2))
  discrap[[region]] <- sdistmatch(a, b)
}

# checking discrepancies
# empty (NULL) list means no discrepancy
lapply(discrap, function(o){ 
  a <- rownames(o)
  b <- o[, 1]
  idx <- a != b
  cbind(a[idx], b[idx])
  })
