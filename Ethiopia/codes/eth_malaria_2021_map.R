
library("sf")
library("tidyverse")
library("readxl")

# =========================================================
# Ethiopian Public Health Institute (EPHI) 
# weekly malaria surveillance data from 2013 to 2022
# restricted access data

# set the path to the malaria data file
data_path <- "~/Downloads/Ethiopia/"

# read the excel file containing the data
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
#eth_data <- eth_data %>%
#  filter(Year <= 2019) 

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
                      "SWEPRS\r\r\n" ~ "South West",
                      .default=RegionName)) %>%
  # change name and spelling of some zones
  # stor modified zone names in the new variable ZoneName2
  mutate(ZoneName2=
           # region by region
           case_match(ZoneName,
                      # Addis Ababa region
                      "Addis Ketema" ~ "Addis Ababa",
                      "Akaki Kaliti" ~ "Addis Ababa",
                      "Arada" ~ "Addis Ababa",
                      "Bole" ~ "Addis Ababa",
                      "Chirkos" ~ "Addis Ababa",
                      "Gulele" ~ "Addis Ababa",
                      "Kolfe Keraniyo" ~ "Addis Ababa",
                      "Lideta" ~ "Addis Ababa",
                      "Nefas Silk Lafto" ~ "Addis Ababa",
                      "Yeka" ~ "Addis Ababa",
                      "Lemi Kura\r\r\n" ~ "Addis Ababa",
                      "Lemi Kura" ~ "Addis Ababa",
                      # Afar region
                      "Zone 01" ~ "Awsi /Zone 1",
                      "Zone 02" ~ "Kilbati /Zone2",
                      "Zone 03" ~ "Gabi /Zone 3",
                      "Zone 04" ~ "Fanti /Zone 4",
                      "Zone 05" ~ "Hari /Zone 5",
                      # Amhara region
                      "Argoba Special Woreda" ~ "South Wello",
                      "Awi" ~ "Awi",
                      "Bahir Dar Liyu Town" ~ "West Gojam",
                      "Dese Town" ~ "South Wello",
                      "East Gojjam" ~ "East Gojam",
                      "Gonder Town" ~ "Central Gondar",
                      "North Gondar" ~ "North Gondar",
                      "North Shewa" ~ "North Shewa (AM)",
                      "North Wollo" ~ "North Wello",
                      "Oromiya" ~ "Oromia",
                      "South Gonder" ~ "South Gondar",
                      "South Wollo" ~ "South Wello",
                      "Wag Himra" ~ "Wag Hamra",
                      "West Gojjam" ~ "West Gojam",
                      "Centeral  Gondar" ~ "Central Gondar",
                      "West Gondar" ~ "West Gondar",
                      # Benishangul-Gumuz region
                      "Assosa" ~ "Assosa",
                      "Kamashi" ~ "Kamashi",
                      "Maokomo Special" ~ "Mao Komo Special",
                      "Metekel" ~ "Metekel",
                      "SHEKA" ~ "Sheka",
                      # Dire Dawa region
                      "Dredewa" ~ "Dire Dawa urban",
                      # Gambela region
                      "Agnuwak" ~ "Agnewak",
                      "Etang Spe." ~ "Itang Special woreda",
                      "Mejenger" ~ "Majang",
                      "Nuwer" ~ "Nuwer",
                      # Harari region
                      "Harari" ~ "Harari",
                      # Oromia region
                      "Adama Special Town" ~ "East Shewa",
                      "Arsi" ~ "Arsi",
                      "Assela Town" ~ "Arsi",
                      "Bale" ~ "Bale",
                      "Bishoftu Town" ~ "East Shewa",
                      "Borena" ~ "Borena",
                      "Burayu Town" ~ "West Shewa",
                      "Dukem Town" ~ "East Shewa",
                      "East Hararge" ~ "East Hararge",
                      "East Shewa" ~ "East Shewa",
                      "East Wellega" ~ "East Wellega",
                      "Finfine Zuria" ~ "Finfine Special",
                      "Gelan Town" ~ "Finfine Special",
                      "Guji" ~ "Guji",
                      "Horo Gudru Wellega" ~ "Horo Gudru Wellega",
                      "Ilu Aba Bora" ~ "Ilu Aba Bora",
                      "Jimma" ~ "Jimma",
                      "Jimma Spe Town" ~ "Jimma",
                      "Lege Dadi Lege Tafo Town" ~ "Finfine Special",
                      "Nekemte Town" ~ "East Wellega",
                      "North Shewa" ~ "North Shewa (OR)",
                      "Qeleme Wellega" ~ "Kelem Wellega",
                      "Sebeta Town" ~ "Finfine Special",
                      "Shashamane Town" ~ "West Arsi",
                      "South West Shewa" ~ "South West Shewa",
                      "Sululta Town" ~ "Finfine Special",
                      "West Arsi" ~ "West Arsi",
                      "West Hararge" ~ "West Hararge",
                      "West Shewa" ~ "West Shewa",
                      "West Wellega" ~ "West Wellega",
                      "Buno Bedele" ~ "Buno Bedele",
                      "west Guji" ~ "West Guji",
                      "Ambo town" ~ "West Shewa",
                      "Holeta town" ~ "Finfine Special",
                      "Modjo town" ~ "East Shewa",
                      "Robe town" ~ "Bale",
                      "Woliso town" ~ "South West Shewa",
                      # Somali region
                      "Afder" ~ "Afder",
                      "Degehabur" ~ "Jarar",
                      "Fik" ~ "Nogob",
                      "Gode" ~ "Shabelle",
                      "Jijiga" ~ "Fafan",
                      "Korahe" ~ "Korahe",
                      "Liben" ~ "Liban",
                      "Shinile" ~ "Siti",
                      "Warder" ~ "Doolo",
                      "Doollo" ~ "Doolo",
                      "FAAFAN" ~ "Fafan",
                      "Jarar" ~ "Jarar",
                      "NOGOB" ~ "Nogob",
                      "SHABEELE" ~ "Shabelle",
                      "SITTI" ~ "Siti",
                      "Dhewa" ~ "Daawa",
                      "Erar" ~ "Erer",
                      "Nogob" ~ "Nogob",
                      # SNNP region
                      "Basketo Town" ~ "Basketo",
                      "Bench Maji" ~ "Bench Sheko",
                      "Dawuro" ~ "Dawuro",
                      "Gamo Gofa" ~ "Gamo",
                      "Gedeo" ~ "Gedeo",
                      "Gurage" ~ "Guraghe",
                      "Hadiya" ~ "Hadiya",
                      "Halaba" ~ "Halaba",
                      "Hawassa Town" ~ "Sidama",
                      "Kefa" ~ "Kefa",
                      "Kembata Tembaro" ~ "Kembata Tembaro",
                      "Konta Town" ~ "Konta Special",
                      "Segen" ~ "Burji",
                      "Sheka" ~ "Sheka",
                      "Sidama" ~ "Sidama",
                      "Siliti" ~ "Siltie",
                      "South Omo" ~ "South Omo",
                      "Wolayita" ~ "Wolayita",
                      "Yem Town" ~ "Yem Special",
                      "Alle" ~ "Alle",
                      "Amaro" ~ "Amaro",
                      "Basketo" ~ "Basketo",
                      "Burji" ~ "Burji",
                      "Derashe" ~ "Derashe",
                      "Gamo" ~ "Gamo",
                      "Konso" ~ "Konso",
                      "Silte" ~ "Siltie",
                      "Yem" ~ "Yem Special",
                      # Sidama region
                      "Hawassa Town" ~ "Sidama",
                      "Hawella Tulla" ~ "Sidama",
                      "Sidama  ZHD" ~ "Sidama",
                      # South West region 
                      "Bench-Sheko" ~ "Bench Sheko",
                      "Dawro  ZHD" ~ "Dawuro",
                      "Kaffa" ~ "Kefa",
                      "KONTA" ~ "Konta Special",
                      "SHEKA" ~ "Sheka",
                      "West OMO ZHD" ~ "Mirab Omo",
                      # Tigray region
                      "Central Tigray" ~ "Central",
                      "Eastern Tigray" ~ "Eastern",
                      "Mekele Especial Zone" ~ "Mekelle",
                      "North Western Tigray" ~ "North Western",
                      "South East" ~ "South Eastern",
                      "South Tigray" ~ "Southern",
                      "Western Tigray" ~ "Western",
                      .default=ZoneName)) %>%
  # fixing renaming of North Shewa zone in two regions
  mutate(ZoneName2=case_when(
    RegionName == "Amhara" & ZoneName == "North Shewa" ~ "North Shewa (AM)",
    RegionName == "Oromia" & ZoneName == "North Shewa" ~ "North Shewa (OR)",
    TRUE ~ ZoneName2
  )) %>%
  # fixing separation of Sidama and South West regions from SNNP
  mutate(RegionName2=case_when(
    RegionName == "SNNP" & ZoneName2 == "Sidama" ~ "Sidama",
    RegionName == "SNNP" & ZoneName2 == "Bench Sheko" ~ "South West",
    RegionName == "SNNP" & ZoneName2 == "Dawuro" ~ "South West",
    RegionName == "SNNP" & ZoneName2 == "Kefa" ~ "South West",
    RegionName == "SNNP" & ZoneName2 == "Konta Special" ~ "South West",
    RegionName == "SNNP" & ZoneName2 == "Sheka" ~ "South West",
    TRUE ~ RegionName
  ))
  # fixing establishment of South West region

# old and new names of regions
eth_data %>% count(RegionName, RegionName2)  %>% print(n=500)

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
for (region in unique(eth_map$ADM1_EN))
{
  cat(region, ": ")
  a <- unique(eth_data %>% 
                filter(RegionName2 == region) %>% 
                pull(ZoneName2))
  b <- unique(eth_map %>% 
                filter(ADM1_EN == region) %>% 
                pull(ADM2_EN))
  discrap[[region]] <- sdistmatch(a, b)
  cat(" done\n")
}

# checking discrepancies
# empty (NULL) list means no discrepancy
lapply(discrap, function(o){ 
  a <- rownames(o)
  b <- o[, 1]
  idx <- a != b
  cbind(a[idx], b[idx])
  })

# which zones on the map are not covered in the data
unique(eth_map %>% 
         pull(ADM2_EN))[ !
                           unique(eth_map %>% pull(ADM2_EN)) %in% unique(eth_data %>% pull(ZoneName2))
         ]

# =========================================================
# exploring missing data

# number of weeks with recorded malaria by zones
eth_data %>% 
  count(Year, Epidemic_Week, RegionName2, ZoneName2) %>%
  count(RegionName2, ZoneName2) %>% arrange(n) %>% print(n=100)

# by region
expand_grid(RegionName2=unique(eth_data %>% pull(RegionName2)),
            Year=2013:2022, Epidemic_Week=1:52) %>%
  left_join(
    eth_data %>% 
      count(RegionName2, Year, Epidemic_Week)
  ) %>%
  filter(is.na(n)) %>%
  select(-n) %>% 
  write_csv(file="~/Desktop/missing_records_region.csv")

# by zone
expand_grid(RegionName2=unique(eth_data %>% pull(RegionName2)),
            ZoneName2=unique(eth_data %>% pull(ZoneName2)),
            Year=2013:2022, Epidemic_Week=1:52) %>%
  left_join(
    eth_data %>% 
      count(RegionName2, ZoneName2, Year, Epidemic_Week)
  ) %>%
  filter(is.na(n)) %>%
  select(-n) %>% 
  write_csv(file="~/Desktop/missing_records_zone.csv")

# =========================================================
# aggregate by zones
eth_data <- eth_data %>% 
  mutate(date2=as.Date(date2)) %>%
  group_by(RegionName2, ZoneName2, date2, Year, Epidemic_Week) %>%
  summarise(
    `Total Malaria Confirmed and Clinical`=
      sum(`Total Malaria Confirmed and Clinical`),
    `TMalaria_OutP_Cases`=
      sum(`TMalaria_OutP_Cases`),
    `TMalaria_InP_Cases`=
      sum(`TMalaria_InP_Cases`),
    `TMalaria_InP_Deaths`=
      sum(`TMalaria_InP_Deaths`),
    `TMSuspected Fever Examined`=
      sum(`TMSuspected Fever Examined`),
    `PosMalaria_RDT_or_Microscopy_PF_OutP_Cases`=
      sum(`PosMalaria_RDT_or_Microscopy_PF_OutP_Cases`),
    `PosMalaria_RDT_or_Microscopy_PV_OutP_Cases`=
      sum(`PosMalaria_RDT_or_Microscopy_PV_OutP_Cases`)
  ) %>%
  rename(RegionName=RegionName2,
         ZoneName=ZoneName2,
         Date=date2)

# save the data as an R data.frame
saveRDS(eth_data, file=paste0(data_path, "eth_data.rds"))
