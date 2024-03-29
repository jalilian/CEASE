
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
if (FALSE)
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
  mutate(An_stephensi_invasive=0)
for (i in 1:nrow(eth_map))
{
  yy <- st_filter(dat, eth_map[i, ]) %>% pull(YEAR_START)
  if (length(yy) > 0)
  {
    idx <- (eth_data$ZoneName == eth_map$ZoneName[i]) &
      (eth_data$Year >= min(yy))
    eth_data$An_stephensi_invasive[idx] <- 1  
  }
}

eth_data %>% group_by(ZoneName) %>% count(An_stephensi_invasive)

rm(dat, spat_covars, spat_temp_covars, i, idx, yy)

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
    xlab('date') + ylab('random effect representing overall long-term trend') +
    #labs(title = NULL, subtitle = NULL, caption = NULL) + 
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
  reffect$idx_zone <- rep(1:92, each=53)
  reffect$ZoneName <- factor(rep(eth_map$ZoneName, each=53))
  ggplot(data=reffect, 
         mapping=aes(x=ID, y=mean)) +
    geom_hline(yintercept=0, linetype="dashed", color="red") +
    geom_line(col="blue") + 
    geom_ribbon(aes(ymin=q0.025, ymax=q0.975), alpha=0.25) +
    xlab('epidemic week') + 
    ylab('random effect representing recurring seasonal patterns') +
    theme_light() + facet_wrap(~ ZoneName, ncol=8)
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
    f(An_stephensi_invasive + 1, An_stephensi_invasive, model="iid") + 
    f(idx_week, model="rw1") +
    f(Epidemic_Week, model="rw2", 
      scale.model=TRUE, constr=TRUE, cyclic=TRUE, 
      group=idx_zone) + #, 
      #control.group=list(model="besag", graph=H, scale.model=TRUE,
      #                   adjust.for.con.comp=TRUE)),
    #f(idx_zone, model="iid"),
    f(idx_zone, model="bym2", graph=H, scale.model=TRUE,
      adjust.for.con.comp=TRUE, constr=TRUE),
  data=eth_data, E=Total_pop,
  family = "nbinomial",
  control.family = list(variant=0),
  control.predictor=list(compute=TRUE, link=1),
  control.compute=list(config=FALSE, waic=TRUE, dic=TRUE, cpo=TRUE,
                       return.marginals.predictor=TRUE),
  num.threads=4, verbose=TRUE, keep=FALSE
  )

printmodel(fit)
pitfun(fit)

fit$summary.hyperpar

trendplot(fit)
seasplot(fit)
spatialplot(fit, what="random_effect")
spatialplot(fit, what="random_effect_exp")
spatialplot(fit, what="linear_predictor")
spatialplot(fit, what="linear_predictor_exp")

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
  inla.zmarginal(inla.tmarginal(exp, a, method = "quantile"))
})

aa <- lapply(fit$marginals.fixed, function(a){
  #inla.emarginal(exp, a)
  x <- exp(inla.rmarginal(10240, a))
  as.list(c(mean(x, trim=0.1), sd(x),
       quantile(x,c(0.025, 0.25, 0.5, 0.75, 0.975))))
})

bb <- round(matrix(unlist(aa), ncol=7, byrow=TRUE)[, c(1, 3, 7)], 3)
rownames(bb) <- names(aa)
bb %>% as.data.frame() %>% rownames_to_column(var="term") %>%
  filter(term != "(Intercept)") %>%
  filter(V2 > 1 | V3 < 1) %>%
  ggplot(aes(x=V1, y=term)) +
  geom_vline(xintercept=1, linetype="dashed", 
             color = "red") +
  geom_point() + 
  geom_errorbar(aes(xmin=V2, xmax=V3)) +
  labs(x="relative risk", y="environmental variables with sifnificant effect") +
  theme_light()



v <- 10000 * fit$summary.fitted.values
v$Year <- eth_data$Year
v$Epidemic_Week <- eth_data$Epidemic_Week
v$idx_week <- eth_data$idx_week
v$ZoneName <- eth_data$ZoneName

g1 <- trendplot(fit)
g2 <- v %>% group_by(Year, Epidemic_Week, idx_week) %>% 
  summarise(mean=mean(mean),
            `0.025quant`=mean(`0.025quant`),
            `0.975quant`=mean(`0.975quant`)) %>%
  mutate(date=as.Date(paste(Year, "01", "01", sep="-")) %m+% 
           weeks(Epidemic_Week) - 1) %>%
  ggplot(aes(x=date, y=mean)) + 
  geom_line(col='blue', linewidth=1.1) + 
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), 
              alpha=0.25) +
  labs(x="date", y="estimated risk of malaria per 10,000 population") +
  theme_light()

ggarrange(g1, g2, ncol=1)


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
  labs(x='epidemic week', y='estimated risk of malaria per 10,000 population') +
  theme_light() + facet_wrap(~ ZoneName, ncol=6)


g1 <- eth_map %>% 
  mutate(idx_zone=1:length(ZoneName)) %>%
  left_join(fit$summary.random$idx_zone %>%
              rename(idx_zone=ID),
            by=join_by(idx_zone)) %>%
  ggplot(aes(fill=mean)) + 
  geom_sf() +
  scale_fill_gradient2(name="spatial random effect", 
                       low='green', 
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
  scale_fill_distiller(name="risk per 10,000 population", 
                       palette = "Reds", direction=1)+
  theme_light()

ggarrange(g1, g2, nrow=1, legend="bottom")


# Bayesian R2
# Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). 
# R-squared for Bayesian regression models. The American Statistician.
nsim <- 1000

yhat <- matrix(NA, nrow=nrow(eth_data), ncol=nsim)
for (i in 1:nrow(eth_data))
{
  yhat[i, ] <- 
    inla.rmarginal(nsim, 
                  fit$marginals.fitted.values[[i]]) * eth_data$Total_pop[i]
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
ggplot(data.frame(R2), aes(x=R2)) + 
  geom_histogram() +
  geom_vline(xintercept = mean(R2), col="blue", 
             linetype="dashed", linewidth=1.25) +
  theme_light()

mean(R2)


mean(apply(yhat, 1, var))
mean(apply(sigma2, 1, mean))
