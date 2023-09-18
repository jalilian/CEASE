# Ethiopia
Ethiopia is a country located in the Horn of Africa and has a diverse and dynamic demographic profile.


## Administrative divisions
The administrative divisions in Ethiopia have considerably changed over time, due to the creation of new woredas, the merger of woredas, and the elevation of some woredas to city administrations.

### Stanford digital repository map
A shapefile containing the boundaries of third-level administrative divisions in Ethiopia in 2015 is freely available on the [Stanford Digital Repository](https://purl.stanford.edu/fx138hn5305) for academic and other non-commercial use. 

This shapefile consists of
- 11 first-level administrative divisions (NAME_1), 9 regional states and 2 chartered cities Addis Ababa and Dire Dawa
- 72 second-level administrative divisions (NAME_2) or zones
- 676 third-level administrative divisions (NAME_3) called woredas

The following command in ```R``` provides this map in [```sf```](https://cran.r-project.org/package=sf) format: 
```R
eth_map <- readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2015_Stanford.rds"))
```
![Administrative divisions of Ethiopia in 2015](/Ethiopia/eth_map_2015.png)

### OCHA map
A shapefile containing the boundaries of administrative divisions in Ethiopia in 2021 is publicly available under the Creative Commons Attribution for Intergovernmental Organisations license by United Nations Office for the Coordination of Humanitarian Affairs (OCHA) country office in Ethiopia on [humdata.com](https://data.humdata.org/dataset/cb58fa1f-687d-4cac-81a7-655ab1efb2d0)

The current version of the shapefile consists of
- 13 first-level administrative divisions (ADM1), 11 regional states and 2 chartered cities Addis Ababa and Dire Dawa
- 92 second-level administrative divisions (ADM2) refered to as zones
- 1082 third-level administrative divisions (ADM3) called woredas

To access map and population data of administrative divisions in Ethiopia in ```R```, simply run the following command:
```R
eth_map <- readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/eth_map.rds"))
```

![Administrative divisions of Ethiopia](/Ethiopia/eth_map.png)

Two woredas are represented in this shapefile because of their contested boundaries: Dawe Serer (ET042199) in Oromia region and Quarsadula (ET050905) in Somali region.


## Population data

According to the [2022 Revision of World Population Prospects (UN WPP)](https://population.un.org/wpp/), the total population of Ethiopia (as of 1 July) has change from 97 million in 2013 to 120 million in 2021.

|Year | Population (thousands) | Median age | Growth rate |
| ---- | ---------- | ----- | ---- | 
2013	|  97 084	| 16.7	| 2.70% |
2014	|  99 747	| 17.0	| 2.71% |
2015	|  102 472	| 17.2	| 2.68% |
2016	|  105 293	| 17.4	| 2.75% |
2017	|  108 198	| 17.7	| 2.70% |
2018	|  111 129	| 17.9	| 2.65% |
2019	|  114 121	| 18.1	| 2.66% |
2020	|  117 191	| 18.3	| 2.65% |
2021	|  120 283	| 18.5	| 2.56% |

The United Nations Office for the Coordination of Humanitarian Affairs (OCHA) country office in Ethiopia also provides gender and age disaggregated estimates of administrative divisions of Ethiopia on [humdata.com](https://data.humdata.org/dataset/cod-ps-eth). 

![Population of Ethiopia](/Ethiopia/eth_pop.png)

Tigray region population estimation is based on 2012 (Ethiopian calendar) population breakdown at kebele level from the Tigray Statistical Authority and Population Sizy by Sex, Area and Density by Region, Zone and Woreda, July 2022 from federal Central Statistics Agency. According to this data, the total estimated population of Ethiopia is around 105 164 thousands. 
