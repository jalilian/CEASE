# Sudan
Sudan, officially the Republic of the Sudan, is a country located in Northeast Africa

## Administrative divisions
Sudan is currently divided into 18 states, and each state is divided into a number of districts.

### OCHA map
A shapefile containing the boundaries of administrative divisions in Sudan in 2021 is publicly available under the Creative Commons Attribution for Intergovernmental Organisations license by United Nations Office for the Coordination of Humanitarian Affairs (OCHA) country office in Sudan on [humdata.com](https://data.humdata.org/dataset/cod-ab-sdn)

The current version of the shapefile consists of
- 19 first-level administrative divisions (ADM1)
  + Abyei PCA, Aj Jazirah, Blue Nile, Central Darfur, East Darfur, Gedaref, Kassala, Khartoum, North Darfur, North Kordofan, Northern, Red Sea, River Nile, Sennar, South Darfur, South Kordofan, West Darfur, West Kordofan, White Nile
- 188 second-level administrative divisions (ADM2) refered to as districts

To access map and population data of administrative divisions in Ethiopia in ```R```, simply run the following command:
```R
sdn_map <- readRDS(url("https://github.com/jalilian/CEASE/raw/main/Sudan/sdn_map.rds"))
```

![Administrative divisions of Sudan](/Sudan/images/sdn_map.png)

Two woredas are represented in this shapefile because of their contested boundaries: Dawe Serer (ET042199) in Oromia region and Quarsadula (ET050905) in Somali region.


## Population data

Population data are provided by [WorldPop Hub](https://hub.worldpop.org/doi/10.5258/SOTON/WP00682). 

