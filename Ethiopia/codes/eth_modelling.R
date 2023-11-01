
library("sf")
library("tidyverse")

# =========================================================
# path to the malaria data file
data_path <- "~/Downloads/Ethiopia/"

# reading the malaria data
eth_data <- readRDS(paste0(data_path, "eth_data.rds"))

# reading spatial covariates 
spat_covars <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/spat_covars.rds"))

# reading spatiotemporal covariates
spat_temp_covars <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/spat_temp_covars.rds"))

# join all covariates
spat_temp_covars <- spat_temp_covars %>% 
  left_join(spat_covars, by=join_by(ADM2_EN)) %>%
  rename(ZoneName=ADM2_EN,
         Year=epi_year, 
         Epidemic_Week=epi_week)

# data from a shapefile provided by OCHA
eth_map <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))
# merge Woreds in each zone
eth_map <- eth_map %>%
  group_by(ADM1_EN, ADM2_EN) %>%
  summarise(geometry=st_union(geometry),
            Total_pop=sum(Total)) %>%
  ungroup() %>%
  rename(RegionName=ADM1_EN, ZoneName=ADM2_EN)
# =========================================================

# remove all records of week 53 of the year 2015
eth_data <- eth_data %>%
  filter(!is.na(Date))

eth_data <- eth_data %>% 
  left_join(tibble(Date=seq(as.Date("2013-01-07"), 
                            as.Date("2022-08-15"), 
                            by=7)) %>%
              mutate(idx_time=1:length(Date)),
            by=join_by(Date)) %>%
  left_join(eth_map %>% 
              mutate(idx_zone=1:length(ZoneName)) %>%
              as_tibble() %>%
              select(RegionName, ZoneName, idx_zone, Total_pop),
            by=join_by(RegionName, ZoneName))

eth_data <- eth_data %>%
  rename(Total_confirmed=`Total Malaria Confirmed and Clinical`)

eth_data <- eth_data %>%
  left_join(spat_temp_covars, 
            by=join_by(ZoneName, Year, Epidemic_Week))
# =========================================================
library("INLA")

###
printmodel <- function(obj)
{
  cat("DIC: ", obj$dic$dic, "\tWAIC: ", obj$waic$waic, 
      "\tfailure: ", sum(obj$cpo$failure > 0, na.rm = TRUE), 
      "\tlogcpo: ", sum(log(obj$cpo$cpo), na.rm = TRUE), "\n")
}

pitfun <- function(fit, J = 10)
{
  u <- (0:J) / J 
  Fbu <- numeric(length(u))
  for (i in 1:length(u))
  {
    v <- fit$cpo$pit - u[i]
    p1 <- (v > 0) & (v <= fit$cpo$cpo)
    p2 <- (v <= 0)
    Fu <- (1 - v / fit$cpo$cpo) * p1+ p2
    Fbu[i] <- mean(Fu, na.rm=TRUE)
  }
  
  H <- hist(fit$cpo$pit, breaks=J)
  H$counts <- diff(Fbu)
  par(mar=c(4, 4, 0.5, 0.5))
  plot(H, xlab="normalized PIT values", ylab="relative frequency", main='')
  abline(h=1/J, lty=2, col=2)
}

trendplot <- function(fit)
{
  reffect <- fit$summary.random$idx_time
  names(reffect) <- c("ID", "mean", "sd", "q0.025", "q0.5", "q0.975", "mode", "kld")
  ggplot(data=reffect, mapping=aes(x=unique(fit$.args$data$Date), y=mean)) +
    geom_line() + 
    geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha=0.25, col='blue') +
    xlab('date') + ylab('RW2 temporal random effect') +
    #labs(title = NULL, subtitle = NULL, caption = NULL) + 
    scale_x_date(date_breaks = "12 month" , date_labels = "%b-%Y") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
}

seasplot <- function(fit)
{
  reffect <- fit$summary.random$Epidemic_Week
  names(reffect) <- c("ID", "mean", "sd", "q0.025", "q0.5", "q0.975", "mode", "kld")
  ggplot(data=reffect, 
         mapping=aes(x=unique(fit$.args$data$Epidemic_Week), y=mean)) +
    geom_line() + 
    geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha=0.25, col='blue') +
    xlab('date') + ylab('RW1 temporal random effect')
}


spatialplot <- function(fit, cmap=eth_map)
{
  dat <- fit$.args$data %>% count(ZoneName, idx_zone) %>%
    left_join(fit$summary.random$idx_zone %>%
              rename(idx_zone=ID),
              by=join_by(idx_zone))
  
  ff <- cmap %>% left_join(dat, by=join_by(ZoneName))
  ggplot(data=ff, aes(fill=mean)) +
    geom_sf()+
    scale_fill_distiller(palette = "Reds", direction=1)
  

}

library("spdep")
tmp <- poly2nb(eth_map, row.names=eth_map$ZoneName)
names(tmp) <- attr(tmp, "region.id")
nb2INLA("map.graph", tmp)
H <- inla.read.graph(filename="map.graph")
image(inla.graph2matrix(H), xlab="", ylab="")

fit <- inla(Total_confirmed ~ 
              u10_mean + u10_min + u10_max + u10_sd +
              v10_mean + v10_min + v10_max + v10_sd + 
              lai_hv_mean + lai_hv_min + lai_hv_max+ lai_hv_sd + 
              lai_lv_mean + lai_lv_min + lai_lv_max + lai_lv_sd +
              skt_mean + skt_min + skt_max + skt_sd +
              #sp_mean + sp_min + sp_max + sp_sd,# +
              tp_mean + tp_min + tp_max + tp_sd + 
              swvl1_mean + swvl1_min + swvl1_max + swvl1_sd +
              pop_density_mean + pop_density_sd + 
              land_cover_mode + land_cover_percent +
              elevation_mean + elevation_sd +
              f(idx_time, model="rw2") +
              f(Epidemic_Week, model="rw2", cyclic=TRUE) +
              f(idx_zone, model='bym2', graph=H, 
                scale.model=TRUE, constr=TRUE, 
                adjust.for.con.comp=TRUE),
              data=eth_data, family = "nbinomial",
            control.predictor=list(compute=TRUE, link=1),
            control.compute=list(config=FALSE, waic=TRUE, dic=TRUE, cpo=TRUE),
            num.threads=4, verbose=TRUE, keep=TRUE)


printmodel(fit)
pitfun(fit)

trendplot(fit)
seasplot(fit)
spatialplot(fit)

fit$summary.fixed
