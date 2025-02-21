
library("sf")
library("tidyverse")
library("ggpubr")
# =========================================================

# data administrative boundary data from a shapefile provided by OCHA
eth_map <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/ETH_Admin_2021_OCHA.rds"))

# merge Woreds in each zone
eth_map <- eth_map %>%
  mutate(ADM2_EN=
           case_match(ADM2_EN,
                      "East Bale" ~ "Bale",
                      "Dire Dawa rural" ~ "Dire Dawa urban",
                      .default=ADM2_EN)) %>%
  group_by(ADM1_EN, ADM2_EN) %>%
  summarise(geometry=st_union(geometry),
            Total_pop=sum(Total)) %>%
  ungroup() %>%
  rename(RegionName=ADM1_EN, ZoneName=ADM2_EN)

# read the EPHI weekly malaria surveillance data
eth_data <- 
  readRDS("~/Downloads/Ethiopia/eth_data.rds") %>%
  # rename the response variable
  rename(Total_confirmed=`Total Malaria Confirmed and Clinical`) %>%
  # arrange records by date
  arrange(Date, RegionName)

# check for ,issing values in the response variable
eth_data %>% filter(is.na(Total_confirmed))
# remove spatio-temporal units with missing values for the response variable
eth_data <- eth_data %>%
  filter(!is.na(Total_confirmed))

# create indices for spatial and temporal units
eth_data <- eth_data %>%
  # create indices for temporal units
  mutate(idx_year=match(Year, sort(unique(Year))),
         idx_week=match(paste0(Year, Epidemic_Week),
                        unique(paste0(Year, Epidemic_Week))),
         idx_week2=idx_week) %>%
  # create region and zone index
  left_join(eth_map %>% 
              mutate(idx_region=match(RegionName, unique(RegionName)),
                     idx_zone=1:length(ZoneName)) %>%
              as_tibble() %>%
              select(RegionName, ZoneName, 
                     idx_region, idx_zone, 
                     Total_pop),
            by=join_by(RegionName, ZoneName)) 


# read spatial covariates 
spat_covars <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/spat_covars.rds"))

spat_covars <- spat_covars %>%
  mutate_if(is.table, as.numeric)

# read spatiotemporal covariates
spat_temp_covars <- 
  readRDS(url("https://github.com/jalilian/CEASE/raw/main/Ethiopia/spat_temp_covars.rds"))

# fix Dire Dawa urban lack of covariates
spat_temp_covars <- 
  bind_rows(spat_temp_covars,
            spat_temp_covars %>%
              filter(ADM2_EN == "Dire Dawa rural") %>%
              mutate(ADM2_EN=case_match(ADM2_EN, 
                                        "Dire Dawa rural" ~ "Dire Dawa urban")))

# merge all covariates
spat_temp_covars <- spat_temp_covars %>% 
  left_join(spat_covars, by=join_by(ADM2_EN)) %>%
  rename(ZoneName=ADM2_EN,
         Year=epi_year, 
         Epidemic_Week=epi_week)

# scale covariates
if (TRUE)
{
  scale2 <- function(x)
  {
    #m <- range(x, na.rm=TRUE)
    #(x - m[1]) / (m[2] - m[1]) 
    (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  }
  spat_temp_covars <- spat_temp_covars %>%
    mutate(across(c(u10_mean:pop_density_sd,
                    elevation_mean:elevation_sd), 
                  scale2))
}

# create date for covariates
spat_temp_covars <-
  spat_temp_covars %>%
  mutate(Date=ymd(paste(Year, "01", "01", sep="-")) + 
           weeks(Epidemic_Week))

# check for NA in covariates
spat_temp_covars %>%
  summarise_all(~ sum(is.na(.))) %>% 
  t()

# =========================================================

eth_data <- eth_data %>%
  # merge covariates
  left_join(spat_temp_covars, 
            by=join_by(ZoneName, Year, Epidemic_Week, Date))

# Repeated weeks 
eth_data %>% 
  count(Year, Epidemic_Week) %>% 
  filter(n > 90)

# Account for population change
# considering 2.7 % annual population growth rate
# according to World Banck: 
#   population of Etiopia in 2016 is 105 million
eth_data <- eth_data %>%
  mutate(Total_pop2 = 
           Total_pop * (1 + 2.7 / 100)^(Year - 2021))

eth_data %>% 
  group_by(Year, Epidemic_Week) %>%
  summarise(sum(Total_pop), sum(Total_pop2), sum(Total_confirmed)) %>%
  print(n=700)


# compute the expected number of cases under the null model
# of spatiotemporal homogeneity
eth_data <- eth_data %>%
  mutate(E=Total_pop2 * 
           sum(Total_confirmed, na.rm=TRUE) / sum(Total_pop2))

# =========================================================
#  read invasive species data
library("readxl")

dat <-
  read_xlsx("~/Downloads/Ethiopia/MTM_INVASIVE_VECTOR_SPECIES_20240322.xlsx",
            sheet="Data") %>% 
  filter(VECTOR_SPECIES_COMPLEX == "An. stephensi") %>%
  select(LONGITUDE, LATITUDE, YEAR_START) %>%
  arrange(YEAR_START) %>%
  st_as_sf(., coords=c("LONGITUDE", "LATITUDE"),
           crs=st_crs(eth_map))

eth_data <- eth_data %>%
  mutate(An_stephensi_invasive=0L)
for (i in 1:nrow(eth_map))
{
  yy <- st_filter(dat, eth_map[i, ]) %>% pull(YEAR_START)
  if (length(yy) > 0)
  {
    idx <- (eth_data$ZoneName == eth_map$ZoneName[i]) &
      (eth_data$Year >= min(yy))
    eth_data$An_stephensi_invasive[idx] <- 1L  
  }
}

eth_data %>% group_by(ZoneName) %>% count(An_stephensi_invasive)

rm(dat, spat_covars, spat_temp_covars, i, idx, yy)

# =========================================================
if (FALSE)
{
  # tables
  smfun <- function(a)
  {
    matrix(c(mean(a), sd(a), min(a), 
             quantile(a, c(0.25, 0.5, 0.75)), max(a)),
           nrow=1)
  }
  
  eth_data %>% group_by(Year) %>%
    summarise(aa=smfun(Total_confirmed)) %>%
    print.data.frame()
  eth_data %>% group_by(RegionName) %>%
    summarise(aa=smfun(Total_confirmed)) %>%
    print.data.frame() 
  # plots
  
  # Fig 1
  # temporal pattern
  g1 <- eth_data %>% 
    group_by(Year, Epidemic_Week) %>% 
    summarise(Total_confirmed=sum(Total_confirmed)) %>% 
    mutate(Date=ymd(paste(Year, "01", "01", sep="-")) + 
             weeks(Epidemic_Week)) %>%
    ggplot(aes(x=Date, y=Total_confirmed)) +
    geom_line() +
    scale_x_date("Date", date_labels = "%b %Y", breaks = "year") +
    scale_y_continuous(labels=function(l) ifelse(l <= 9999, l, scales::comma(l))) +
    #stat_smooth(se=FALSE) +
    labs(y="Total number of clinical malaria cases", x="Date") +
    theme_light()
  
  
  library("ggspatial")
  # Fig 2
  # spatial pattern
  g2 <- eth_map %>%
    left_join(eth_data %>% 
                group_by(ZoneName, Year) %>%
                summarise(Total_confirmed=sum(Total_confirmed)) %>%
                pivot_wider(names_from=Year, 
                            values_from=Total_confirmed),
              by=join_by(ZoneName)) %>% 
    pivot_longer(cols=`2013`:`2023`, 
                 names_to="Year", 
                 values_to ="Total_confirmed") %>%
    # account for population changes by considering constant 2.67% population growth rate
    mutate(Total_pop2=
             Total_pop * ( (1 - 2.67 / 100)^(2022 - as.numeric(Year)) )) %>%
    mutate(rate=Total_confirmed / Total_pop2 * 1e4) %>%
    ggplot(aes(fill=rate)) + 
    geom_sf() +
    scale_fill_distiller(name="Rate per 10,000 population", 
                         palette = "Reds", direction=1)+
    facet_wrap(~Year, ncol=4) +
    annotation_north_arrow(location = "tr", 
                           height = unit(0.8, "cm"),
                           width = unit(0.8, "cm")) +
    theme_light() +
    theme(legend.position='bottom',
          legend.key.width = unit(2, 'cm'),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ggarrange(g1, g2, ncol=1, heights=c(1, 2.5))
  ggsave(filename="Fig1_new.pdf", width=8, height=9)
  library("devEMF")
  ggsave(filename="Fig1_new.emf", device=emf, width=9, height=10)
}

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
  reffect <- fit$summary.random$idx_week
  names(reffect) <- c("ID", "mean", "sd", "q0.025", "q0.5", "q0.975", "mode", "kld")
  ggplot(data=reffect, 
         mapping=aes(x=fit$.args$data$Date[1] %m+% 
                       weeks(unique(fit$.args$data$idx_week) - 1),
                     y=mean)) +
    geom_hline(yintercept=0, linetype="dashed", 
               color = "red") +
    geom_line(col='blue', linewidth=1.1) + 
    geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha=0.25) +
    labs(x="", 
         y=expression(random~effect~tau[t]~representing~overall~long-term~trend)) +
    #scale_x_date(date_breaks = "12 month" , date_labels = "%b-%Y") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    theme_light()
}

seasplot <- function(fit)
{
  reffect <- fit$summary.random$Epidemic_Week
  names(reffect) <- c("ID", "mean", "sd", 
                      "q0.025", "q0.5", "q0.975", 
                      "mode", "kld")
  reffect$idx_zone <- rep(1:length(eth_map$ZoneName), each=53)
  reffect$ZoneName <- factor(rep(eth_map$ZoneName, each=53))
  ggplot(data=reffect, 
         mapping=aes(x=ID, y=mean)) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    geom_line(col="blue") + 
    geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha=0.25) +
    xlab('epidemic week') + 
    ylab(expression(random~effect~s[it]~representing~recurring~seasonal~patterns)) +
    theme_light() + facet_wrap(~ ZoneName, ncol=9)
}

seasplot2 <- function(fit)
{
  reffect <- fit$summary.random$Epidemic_Week
  names(reffect) <- c("ID", "mean", "sd", 
                      "q0.025", "q0.5", "q0.975", 
                      "mode", "kld")
  ggplot(data=reffect, 
         mapping=aes(x=unique(fit$.args$data$Epidemic_Week), y=mean)) +
    geom_hline(yintercept=0, linetype="dashed", 
               color = "red") +
    geom_line() + 
    geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha=0.25, col='blue') +
    xlab('Epidemic week') + ylab('Recurring seasonal patterns') +
    theme_light()
}


spatialplot <- function(fit, what="random_effect", cmap=eth_map)
{
  switch (what, random_effect={
    dat <- fit$.args$data %>% count(ZoneName, idx_zone) %>%
      left_join(fit$summary.random$idx_zone %>%
                  rename(idx_zone=ID),
                by=join_by(idx_zone))
    
    midpoint <- 0
  }, random_effect_exp={
    out <- parallel::mclapply(
      fit$marginals.random$idx_zone, 
      function(o){ 
        fun <- function(v)
        {
          exp(v) / (1 + exp(v))
        }
        inla.zmarginal(inla.tmarginal(fun, o),
                       silent=TRUE) 
      }, mc.cores=4)
    
    out <- bind_rows(out) %>%
      bind_cols(idx_zone=fit$summary.random$idx_zone$ID) %>%
      rename("0.025quant"="quant0.025",
             "0.975quant"="quant0.975")
    
    dat <- fit$.args$data %>% count(ZoneName, idx_zone) %>%
      left_join(out,
                by=join_by(idx_zone))
    
    midpoint <- mean(dat$mean)
  }, linear_predictor={
    dat <- bind_cols(fit$.args$data, 
              fit$summary.linear.predictor) %>%
      group_by(RegionName, ZoneName, 
               idx_region, idx_zone) %>%
      summarise(mean=median(mean),
                `0.025quant`=median(`0.025quant`),
                `0.975quant`=median(`0.975quant`)) %>%
      ungroup()
    
    midpoint <- 0
  }, linear_predictor_exp={
    out <- parallel::mclapply(
      fit$marginals.linear.predictor, 
      function(o){ 
        fun <- function(v)
        {
          exp(v) / (1 + exp(v))
        }
        inla.zmarginal(inla.tmarginal(fun, o),
                       silent=TRUE) 
      }, mc.cores=4)
    
    out <- bind_rows(out) %>%
      rename("0.025quant"="quant0.025",
             "0.975quant"="quant0.975")
    
    dat <- bind_cols(fit$.args$data, 
              out) %>%
      group_by(RegionName, ZoneName, 
               idx_region, idx_zone) %>%
      summarise(mean=median(mean),
                `0.025quant`=median(`0.025quant`),
                `0.975quant`=median(`0.975quant`)) %>%
      ungroup()
    
    midpoint <- mean(dat$mean)
    })
  
  ff <- cmap %>% 
    left_join(dat, by=join_by(ZoneName)) %>%
    mutate(label=case_when(
      `0.025quant` > 0 ~ "+",
      `0.975quant` < 0 ~ "-",
      .default = NA
    ))
  
  g0 <- ggplot() +
    geom_sf(data=ff) +
    geom_sf(data=ff %>% filter(`0.025quant` > 0),
            mapping=aes(fill=mean)) + 
    geom_sf(data=ff %>% filter(`0.975quant` < 0),
            mapping=aes(fill=mean)) +
    geom_sf(data=ff %>% filter(`0.975quant` > 0 &
                                 `0.025quant` < 0),
            mapping=aes(fill=0)) +
    scale_fill_gradient2(low='green', 
                         mid='white', 
                         high='red', 
                         midpoint = midpoint,
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "fill") +
    theme_light()
  
  g1 <- ggplot() +
    geom_sf(data=ff, mapping=aes(fill=`0.025quant`)) +
    scale_fill_distiller(palette = "Greens", direction=1) +
    theme_light()
  g2 <- ggplot() +
    geom_sf(data=ff, mapping=aes(fill=`0.975quant`)) +
    scale_fill_distiller(palette = "Reds", direction=1) +
    theme_light()
  ggarrange(g0, ggarrange(g1, g2, nrow=1), ncol=1)
}

library("spdep")
tmp <- poly2nb(eth_map, row.names=eth_map$ZoneName)
names(tmp) <- attr(tmp, "region.id")
nb2INLA("map.graph", tmp)
H <- inla.read.graph(filename="map.graph")
image(inla.graph2matrix(H), xlab="", ylab="")


fit <- inla(
  Total_confirmed ~ 
    # fixed covariate effects
    u10_mean + u10_sd + u10_min + u10_max + 
    v10_mean + v10_sd + v10_min + v10_max + 
    lai_hv_mean +  lai_hv_sd + lai_hv_min + lai_hv_max +
    lai_lv_mean + lai_lv_sd + lai_lv_min + lai_lv_max +
    skt_mean + skt_sd + skt_min + skt_max +
    #sp_mean + sp_sd + #sp_min + sp_max + 
    tp_mean + tp_sd + tp_min + tp_max + 
    swvl1_mean + swvl1_sd + swvl1_min + swvl1_max + 
    pop_density_mean + pop_density_sd + pop_density_min + pop_density_max +
    land_cover_1 + land_cover_2 + land_cover_3 +
    land_cover_4 + land_cover_5 + land_cover_6 +
    land_cover_7 + land_cover_8 + land_cover_10 +
    elevation_mean + elevation_sd + elevation_min + elevation_max +
    # random spatial and temporal effects
    #f(An_stephensi_invasive + 1, An_stephensi_invasive, model="iid") +
    f(An_stephensi_invasive, model="iid") +
    f(idx_week, model="rw1") +
    f(Epidemic_Week, model="rw2", scale.model=TRUE, constr=TRUE, 
      cyclic=TRUE, group=idx_zone) + #, 
      #control.group=list(model="besag", graph=H, scale.model=TRUE,
      #                   adjust.for.con.comp=TRUE)),
   # f(idx_zone, model="iid"),
    f(idx_zone, model="bym", 
      graph=H, scale.model=TRUE, adjust.for.con.comp=TRUE, 
      constr=TRUE),
  data=eth_data, E=Total_pop2,
  family = "nbinomial",
  control.family = list(variant=0),
  control.fixed=list(mean=0, prec=1e-06, correlation.matrix=TRUE),
  control.predictor=list(compute=TRUE, link=1),
  control.compute=list(cofig=FALSE, waic=TRUE, dic=TRUE, cpo=TRUE,
                       return.marginals.predictor=TRUE),
  control.inla = list(strategy = "laplace", npoints = 27, 
                      int.strategy = "ccd", diff.logdens = 3),
  num.threads=2, verbose=TRUE, keep=FALSE
  )

printmodel(fit)
pitfun(fit)

saveRDS(fit, file="~/Downloads/Ethiopia/eth_fit.rds")

fit$summary.hyperpar

# variance components
g1 <- do.call(rbind.data.frame, 
        lapply(fit$marginals.hyperpar[c(2:4, 6:7)], 
               function(v){
                 c(inla.zmarginal(inla.tmarginal(
                   function(o){ -log( (o + 0.5 +.Machine$double.eps)) }, v,
                   n = 10 * 2048L,),
                   silent=TRUE)[c(1, 3, 7)])
                 })) %>%
  mutate(component=rownames(.)) %>%
  mutate(component=
           case_match(component,
                      "Precision for An_stephensi_invasive" ~ 
                        "\u03b8: invasive species",
                      "Precision for idx_week" ~ 
                        "\u03C4: temporal trend",
                      "Precision for Epidemic_Week" ~ 
                        "s: zone-specific seasonality",
                      "Precision for idx_zone (iid component)" ~ 
                        "\u03be: spatial (iid component)",
                      "Precision for idx_zone (spatial component)" ~ 
                        "\u03be: spatial (structured component)")) %>%
  as_tibble() %>%
  ggplot(aes(y=component, x=mean)) +
  geom_point() +
  geom_errorbar(aes(xmin=`quant0.025`, xmax=`quant0.975`)) +
  labs(y="Considered random effects in the model", x="logarithm of variance") +
  theme_light() 


fit$summary.hyperpar[c(2:4, 6:7), ] %>%
  mutate(component=rownames(.)) %>%
  mutate(component=
           case_match(component,
                      "Precision for An_stephensi_invasive" ~ 
                        "Precision for \u0398",
                      "Precision for idx_week" ~ 
                        "Precision for \u03c4",
                      "Precision for Epidemic_Week" ~ 
                        "Precision for s",
                      "Precision for idx_zone (iid component)" ~ 
                        "Precision for \u03be (iid component)",
                      "Precision for idx_zone (spatial component)" ~ 
                        "Precision for \u03be (spatial component)")) %>%
  ggplot(aes(y=component, x=mean)) +
  geom_point() +
  geom_errorbar(aes(xmin=`0.025quant`, xmax=`0.975quant`)) +
  labs(y="precision of random effects", x="") +
  theme_light()


trendplot(fit)
seasplot(fit)
spatialplot(fit, what="random_effect")
spatialplot(fit, what="random_effect_exp")
spatialplot(fit, what="linear_predictor")
spatialplot(fit, what="linear_predictor_exp")
library("devEMF")
ggsave(filename="Fig4_supp.emf", device=emf, width=10, height=8)


fit$summary.fixed %>%
  rownames_to_column(var="term") %>%
  filter(term != "(Intercept)") %>%
  filter(`0.025quant` * `0.975quant` > 0) %>% 
  select(term, mean, `0.025quant`, `0.975quant`) %>%
  ggplot(aes(x=mean, y=term)) +
  geom_vline(xintercept=0, 
             linetype="dashed", 
             color = "red") +
  geom_point() + 
  geom_errorbar(aes(xmin=`0.025quant`, xmax=`0.975quant`)) +
  theme_light()


aa <- lapply(fit$marginals.fixed, function(a){
  #inla.emarginal(exp, a)
  x <- exp(inla.rmarginal(10240, a))
  as.list(c(mean(x, trim=0.1), sd(x),
       quantile(x,c(0.025, 0.25, 0.5, 0.75, 0.975))))
})

bb <- round(matrix(unlist(aa), ncol=7, byrow=TRUE)[, c(1, 3, 7)], 3)
rownames(bb) <- names(aa)
g2 <- bb %>% as.data.frame() %>% rownames_to_column(var="term") %>%
  filter(term != "(Intercept)") %>%
  filter(V2 > 1 | V3 < 1) %>%
  ggplot(aes(x=V1, y=term)) +
  geom_vline(xintercept=1, linetype="dashed", 
             color = "red") +
  geom_point() + 
  geom_errorbar(aes(xmin=V2, xmax=V3)) +
  labs(x="relative risk", y="Environmental variables with sifnificant effect") +
  theme_light() + 
  scale_y_discrete(  position = "right") 

library("ggpubr")
library("devEMF")
emf(file="~/Fig6.emf", width=9, height=4.5)
ggarrange(g1, g2, labels=c("A", "B"),  nrow=1, ncol=2)
dev.off()

# Fig 6
tiff("~/Fig6.tiff", width=3 * 950, height=3 * 400, res=300)
ggarrange(g1, g2, labels=c("A", "B"),  nrow=1, ncol=2)
dev.off()

v <- 10000 * fit$summary.fitted.values
v$Year <- eth_data$Year
v$Epidemic_Week <- eth_data$Epidemic_Week
v$idx_week <- eth_data$idx_week
v$ZoneName <- eth_data$ZoneName

g1 <- trendplot(fit)
g2 <- v %>% group_by(Year, Epidemic_Week, idx_week) %>% 
  summarise(mean=sum(mean),
            `0.025quant`=sum(`0.025quant`),
            `0.975quant`=sum(`0.975quant`)) %>%
  mutate(date=as.Date(paste(Year, "01", "01", sep="-")) %m+% 
           weeks(Epidemic_Week) - 1) %>%
  ggplot(aes(x=date, y=mean)) + 
  geom_line(col='blue', linewidth=1.05) + 
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), 
              alpha=0.25) +
  labs(x="", y=expression(aggregated~risk~theta[t]~per~10000~population)) +
  theme_light()

ggarrange(g1, g2, ncol=1)

# Fig 3
tiff("~/Fig3.tiff", width=3 * 950, height=3 * 400, res=300)
g2
dev.off()

v %>% group_by(Epidemic_Week) %>% summarise(a=mean(mean)) %>%
  ggplot(aes(x=Epidemic_Week, y=a)) + geom_line()

v %>% group_by(Epidemic_Week, ZoneName) %>%
  summarise(mean=mean(mean), 
            `0.025quant`=mean(`0.025quant`),
            `0.975quant`=mean(`0.975quant`)) %>%
  ggplot(mapping=aes(x=Epidemic_Week, y=mean)) +
  #geom_hline(yintercept=0, linetype="dashed", color="red") +
  geom_line(col='blue') + 
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), 
              alpha=0.25) +
  labs(x='epidemic week', 
       y=expression(risk~of~malaria~theta[it]~per~10000~population)) +
  theme_light() + facet_wrap(~ ZoneName, ncol=9)


g1 <- eth_map %>% 
  mutate(idx_zone=1:length(ZoneName)) %>%
  left_join(fit$summary.random$idx_zone %>%
              rename(idx_zone=ID),
            by=join_by(idx_zone)) %>%
  ggplot(aes(fill=mean)) + 
  geom_sf() +
  annotation_north_arrow(location = "tr") +
  scale_fill_gradient2(name="Magnitude and sign of the spatial random effect", 
                       low='blue', 
                       mid='white', 
                       high='red', 
                       midpoint = 0,
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "fill") +
  theme_light()

g2 <- eth_map %>% left_join(
  v %>% group_by(ZoneName) %>% 
    summarise(mean=mean(mean)),
  by = join_by(ZoneName)
) %>%
  ggplot(aes(fill=mean)) + 
  geom_sf() +
  annotation_north_arrow(location = "tr") +
  scale_fill_distiller(name="Malaria risk per 10,000 population", 
                       palette = "Reds", direction=1)+
  theme_light()

# Fig 5
emf(file="~/Fig5.emf", width=9, height=4.15)
ggarrange(g1, g2, nrow=1, legend="bottom")
dev.off()

tiff("~/Fig5.tiff", width=3 * 950, height=3 * 450, res=300)
ggarrange(g1, g2, nrow=1, legend="bottom")
dev.off()

# cluster analysis
zn <- eth_map$ZoneName
u <- fit$summary.random$Epidemic_Week[, 2:7]
u$ZoneName <- factor(rep(zn, each=53))
dd <- matrix(0, nrow=length(zn), ncol=length(zn))
for (i in 2:90)
{
  for (j in 1:(i- 1))
  {
    dd[i, j] <- sum((u[u$ZoneName == zn[i], -7] - u[u$ZoneName == zn[j], -7])^2)
    dd[j, i] <- dd[i, j]
  }
}
colnames(dd) <- rownames(dd) <- zn
hc <- hclust(as.dist(dd))
plot(hc)

hcc <- cutree(hc, k=4)

library("ggspatial")
# Fig 4
tiff("~/Fig4.tiff", width=3 * 975, height=3 * 675, res=300)
g0 <- eth_map %>% 
  left_join(data.frame(ZoneName=names(hcc), 
                       cluster=factor(paste("cluster", unname(hcc))))) %>%
  ggplot(aes(fill=cluster)) + 
  geom_sf() +
  annotation_north_arrow(location = "tr") +
  #geom_sf_text(aes(label=ZoneName), 
  #             position = position_jitter(width=0.2, height=0.2), 
  #             size=2, angle=0, fontface = "bold") +
  labs(x="", y="", fill="Zone clusters") +
  theme_light()
g0
dev.off()

# clusters by risk
eth_data %>% 
  mutate(risk=fit$summary.fitted.values$mean * 1e4,
         risk0.025=fit$summary.fitted.values$`0.025quant` * 1e4,
         risk0.975=fit$summary.fitted.values$`0.975quant` * 1e4) %>%
  group_by(ZoneName) %>% 
  summarise(risk=mean(risk), 
            risk0.025=mean(risk0.025), 
            risk0.975=mean(risk0.975)) %>%
  ungroup() %>%
  left_join(data.frame(ZoneName=names(hcc), cluster=factor(unname(hcc)))) %>%
  group_by(cluster) %>% 
  summarise(mean(risk), mean(risk0.025), mean(risk0.975))


g1 <- u %>% 
  mutate(Epidemic_Week=rep(1:53, 90)) %>%
  filter(ZoneName %in% names(hcc[hcc == 1])) %>%
  ggplot(aes(x=Epidemic_Week, y=mean, group=ZoneName)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`),
              alpha=0.2, col=NA, fill="blue") +
  geom_line(col="navy", alpha=0.5) +
  geom_hline(aes(yintercept = 0), col="red", linetype="dashed") +
  labs(x="epidemic week", y="seasonal effect") +
  theme_light()

g2 <-  u %>% 
  mutate(Epidemic_Week=rep(1:53, 90)) %>%
  filter(ZoneName %in% names(hcc[hcc == 2])) %>%
  ggplot(aes(x=Epidemic_Week, y=mean, group=ZoneName)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`),
              alpha=0.2, col=NA, fill="blue") +
  geom_line(col="navy", alpha=0.5) +
  geom_hline(aes(yintercept = 0), col="red", linetype="dashed") +
  labs(x="epidemic week", y="seasonal effect") +
  theme_light()

g3 <-  u %>% 
  mutate(Epidemic_Week=rep(1:53, 90)) %>%
  filter(ZoneName %in% names(hcc[hcc == 3])) %>%
  ggplot(aes(x=Epidemic_Week, y=mean, group=ZoneName)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`),
              alpha=0.2, col=NA, fill="blue") +
  geom_line(col="navy", alpha=0.5) +
  geom_hline(aes(yintercept = 0), col="red", linetype="dashed") +
  labs(x="epidemic week", y="seasonal effect") +
  theme_light()

g4 <-  u %>% 
  mutate(Epidemic_Week=rep(1:53, 90)) %>%
  filter(ZoneName %in% names(hcc[hcc == 4])) %>%
  ggplot(aes(x=Epidemic_Week, y=mean, group=ZoneName)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`),
              alpha=0.2, col=NA, fill="blue") +
  geom_line(col="navy", alpha=0.5) +
  geom_hline(aes(yintercept = 0), col="red", linetype="dashed") +
  labs(x="epidemic week", y="seasonal effect") +
  theme_light()

ggarrange(g1, g2, g3, g4, ncol=2, nrow=2,
          labels=paste("Cluster", 1:4), label.x=0.05)

emf(file="~/Fig4.emf", width=7, height=8.5)
ggarrange(g0, 
  ggarrange(g1 + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank()), 
            g2 + theme(axis.title.y=element_blank()) +
              theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank()), 
            g3, 
            g4 + theme(axis.title.y=element_blank()), 
            ncol=2, nrow=2,
            labels=paste("Cluster", 1:4), label.x=0.05), 
  ncol=1, nrow=2)
dev.off()

u %>%
  left_join(data.frame(hcc, ZoneName=names(hcc)),
            by=join_by(ZoneName)) %>% 
  mutate(Epidemic_Week=rep(1:53, 90)) %>%
  ggplot(aes(x=Epidemic_Week, y=mean, group=ZoneName)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`),
              alpha=0.2, col=NA, fill="blue") +
  geom_line(col="navy", alpha=0.5) +
  geom_hline(aes(yintercept = 0), col="red", linetype="dashed") +
  labs(x="epidemic week", y="seasonal effect") +
  theme_light() + 
  facet_wrap(~hcc, ncol=1) + 
  theme(strip.background =element_rect(fill="blue4"))

#library(ggdendro)
#ggdendrogram(hc, rotate = FALSE, size = 2)
# Bayesian R2
# Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). 
# R-squared for Bayesian regression models. The American Statistician.

R2fun <- function(fit, fit0=NULL, nsim=5000)
{
  n <- nrow(fit$.args$data)
  E <- fit$.args$E
  yhat <- matrix(NA, nrow=n, ncol=nsim)
  for (i in 1:n)
  {
    yhat[i, ] <- 
      inla.rmarginal(nsim, fit$marginals.fitted.values[[i]]) * E[i]
  }
  # explained variance
  varfit <- apply(yhat, 2, var)
  size <- 
    inla.rmarginal(nsim, 
                   fit$marginals.hyperpar$`size for the nbinomial observations (1/overdispersion)`)
  # size <- rep(fit$summary.hyperpar$mean[1], nsim)
  # residual variance
  sigma2 <- matrix(NA, nrow=nrow(eth_data), ncol=nsim)
  for (j in 1:nsim)
  {
    # negative binomial model
    sigma2[, j] <- yhat[, j] * (1 + yhat[, j] / size[j])
  }
  varres <- apply(sigma2, 2, mean)
  
  R2 <- varfit / (varfit + varres)
  if (!is.null(fit0))
  {
    yhat0 <- matrix(NA, nrow=n, ncol=nsim)
    for (i in 1:n)
    {
      yhat0[i, ] <- 
        inla.rmarginal(nsim, fit0$marginals.fitted.values[[i]]) * E[i]
    }
    varrand <- apply(yhat - yhat0, 2, var)
    R2partial <- varrand / (varrand + varres)
    attr(R2, "R2partial") <- R2partial
  }
  return(R2)
}

R2 <- R2fun(fit, fit0)
ggplot(data.frame(R2), aes(x=R2)) + 
  geom_histogram(aes(y = after_stat(count / sum(count)))) +
  geom_vline(xintercept = mean(R2), col="blue", 
             linetype="dashed", linewidth=1.25) +
  labs(x=expression(R^{2}), y="probability density") + 
  theme_light()

mean(R2)
summary(attr(R2, "R2partial"))


library("terra")
library("geodata")
library("tidyterra")
library("ggspatial")
library("devEMF")
library("ggpubr")

eth_alt <- elevation_30s(country="ETH", path=tempdir())
eth_alt_cat <- classify(eth_alt, c(-Inf, 1500, 2400, Inf))
levels(eth_alt_cat) <- data.frame(ID=0:2, 
                                  Altitude=c("Kola", "Weyna Dega", "Dega"))

ggplot() +
  geom_spatraster(data=eth_alt_cat, na.rm=FALSE) + 
  scale_fill_discrete(na.translate=FALSE) +
  geom_sf(data=eth_map, fill=NA) +
  annotation_north_arrow(location = "tr") +
  theme_light()
ggsave(filename="~/ETH_altitude.emf", device=emf, width=10, height=8)
