
setwd("~/Downloads/Ethiopia/ITS/")

library("tidyverse")

# read the data
mata <- read_csv("Metehara.csv") %>%
  mutate(Town=case_match(Town,
                         "Metehara Town" ~ "Matahara",
                         .default=Town))
dire <- read_csv("Dire Dawa.csv")
batu <- read_csv("Batu.csv")

# imputation of missing data using moving averages
zerofun <- function(dat, k=3)
{
  y <- dat$RR # dat$`Confirmed malaria cases`
  y[y == 0] <- NA
  wt <- c(1/(k:1), 1/(1:k))
  for (i in 1:nrow(dat))
  {
    if (is.na(y[i]))
    {
      idx <- c(i - 1:k, i + 1:k)
      ok <- (idx >= 1) & (idx <= nrow(dat))
      y[i] <-  weighted.mean(y[idx[ok]], wt[ok], na.rm=TRUE)#mean(y[c(i - 1:k, i + 1:k)], na.rm=TRUE)
      
    }
  }
  dat$rate <- y
  dat$date <- ymd(paste0(dat$Year, "-01-01")) + 
    weeks(dat$Epidemic_Week - 1)
  return(dat)
}

mata <- zerofun(mata)
dire <- zerofun(dire)
batu <- zerofun(batu)

library("ggpubr")
g1 <- batu %>% 
  ggplot(aes(x=date, y=(rate))) +
  geom_line() + 
  geom_vline(xintercept= as.Date("2022-08-13"), 
             col="blue", linetype=2) +
  labs(x="Date", y="Malaria incidence rate per 10,000") +
  theme_classic()
g2 <- dire %>% 
  ggplot(aes(x=date, y=(rate))) +
  geom_line() + 
  geom_vline(xintercept= as.Date("2022-08-13"), 
             col="blue", linetype=2) +
  labs(x="Date", y="Malaria incidence rate per 10,000") +
  theme_classic()
g3 <- mata %>% 
  ggplot(aes(x=date, y=(rate))) +
  geom_line() +
  labs(x="Date", y="Malaria incidence rate per 10,000") +
  theme_classic()

ggarrange(g1, g2, g3, nrow=3, ncol=1,
          labels=c("Batu", "Dire Dawa", "Metehara"),
          label.x=0.05)

library("devEMF")
ggsave(filename="Fig1.emf", device=emf, width=10, height=8)


# modelling the data using ITS
d1 <- data.frame(date=batu$date, y=batu$rate, t=1:nrow(batu)) %>%
  mutate(z=as.numeric(date >= as.Date("2022-08-13")), 
         s=pmax(0, t - max(which(z == 0))),
         tidx=t,
         sidx=batu$Epidemic_Week)
d2 <- data.frame(date=dire$date, y=dire$rate, t=1:nrow(dire)) %>%
  mutate(z=as.numeric(date >= as.Date("2022-08-13")), 
         s=pmax(0, t - max(which(z == 0))),
         tidx=t,
         sidx=dire$Epidemic_Week)
d3 <- data.frame(date=mata$date, y=mata$rate, t=1:nrow(mata))%>%
  mutate(z=as.numeric(date >= as.Date("2022-08-13")), 
         s=pmax(0, t - max(which(z == 0))),
         tidx=t,
         sidx=mata$Epidemic_Week)


library("trend")
# non-parametric Cox and Stuart trend test
# tests trend change in the first third of observations with the last third of observations
cs.test(log(d1$y))
cs.test(log(d2$y))
cs.test(log(d3$y))

# change point detection
library("changepoint")
cpt.mean(log(d1$y), method="PELT")
d1$date[183]
cpt.mean(log(d2$y), method="PELT")
d2$date[c(28, 88, 139)]
cpt.mean(log(d3$y), method="PELT")
d3$date[134]

library("strucchange")
bp1 <- breakpoints(log(y) ~ t, data=d1, h=0.2)
d1$date[bp1$breakpoints]
bp2 <- breakpoints(log(y) ~ t, data=d2, h=0.2)
d2$date[bp2$breakpoints]
bp3 <- breakpoints(log(y) ~ t, data=d3, h=0.2)
d3$date[bp3$breakpoints]


library("INLA")

fm <- log(y) ~ t + z + s + 
  #f(tidx, model="ar1") +
  f(sidx, model="rw2", cyclic=TRUE, constr=TRUE)
fit1 <- inla(fm, data=d1)
round(fit1$summary.fixed[, c(1, 3, 5)], 3)

fit2 <- inla(fm, data=d2)
round(fit2$summary.fixed[, c(1, 3, 5)], 3)

fit3 <- inla(fm, data=d3)
round(fit3$summary.fixed[, c(1, 3, 5)], 3)


fit1$summary.fitted.values %>%
  ggplot(aes(x=d1$date, y=mean)) + 
  geom_line() + 
  geom_line(aes(x=d1$date, y=log(d1$y)))+
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), 
              fill="blue", alpha=0.5) +
  theme_classic()

fit2$summary.fitted.values %>%
  ggplot(aes(x=d2$date, y=mean)) + 
  geom_line() + 
  geom_line(aes(x=d2$date, y=log(d2$y)))+
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), 
              fill="blue", alpha=0.5) +
  theme_classic()

fit3$summary.fitted.values %>%
  ggplot(aes(x=d3$date, y=mean)) + 
  geom_line() + 
  geom_line(aes(x=d3$date, y=log(d3$y)))+
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), 
              fill="blue", alpha=0.5) +
  theme_classic()

mu1 <- as.matrix(cbind(1, d1$t)) %*% 
  as.matrix(fit1$summary.fixed[1:2, ])

nu1 <- as.matrix(cbind(d1$z, d1$s)) %*% 
  as.matrix(fit1$summary.fixed[3:4, ])
nu1[nu1 == 0] <- NA
nu1 <- mu1 + nu1
mu1[!is.na(nu1)] <- NA
  
g1 <- d1 %>% 
  ggplot(aes(x=date, y=log(y))) + 
  geom_line() +
  geom_line(aes(x=date, y=fit1$summary.fitted.values$mean), col="green4", linetype=2) +
  geom_ribbon(aes(ymin=fit1$summary.fitted.values$`0.025quant`,
                  ymax=fit1$summary.fitted.values$`0.975quant`),
              col="green", fill="green", alpha=0.5) + 
  geom_line(aes(x=date, y=mu1[, 1]), col="blue") +
  geom_ribbon(aes(ymin=mu1[, 3], ymax=mu1[, 5]), col="blue", fill="blue", alpha=0.5) +
  geom_line(aes(x=date, y=nu1[, 1]), col="red") +
  geom_ribbon(aes(ymin=nu1[, 3], ymax=nu1[, 5]), col="red", fill="red", alpha=0.5) +
  labs(x="Date", y="Logarithm of malaria incidence rate per 10,000") +
  theme_classic()


mu2 <- as.matrix(cbind(1, d2$t)) %*% 
  as.matrix(fit2$summary.fixed[1:2, ])

nu2 <- as.matrix(cbind(d2$z, d2$s)) %*% 
  as.matrix(fit2$summary.fixed[3:4, ])
nu2[nu2 == 0] <- NA
nu2 <- mu2 + nu2
mu2[!is.na(nu2)] <- NA

g2 <- d2 %>% 
  ggplot(aes(x=date, y=log(y))) + 
  geom_line() +
  geom_line(aes(x=date, y=fit2$summary.fitted.values$mean), col="green4", linetype=2) +
  geom_ribbon(aes(ymin=fit2$summary.fitted.values$`0.025quant`,
                  ymax=fit2$summary.fitted.values$`0.975quant`),
              col="green", fill="green", alpha=0.5) + 
  geom_line(aes(x=date, y=mu2[, 1]), col="blue") +
  geom_ribbon(aes(ymin=mu2[, 3], ymax=mu2[, 5]), col="blue", fill="blue", alpha=0.5) +
  geom_line(aes(x=date, y=nu2[, 1]), col="red") +
  geom_ribbon(aes(ymin=nu2[, 3], ymax=nu2[, 5]), col="red", fill="red", alpha=0.5) +
  labs(x="Date", y="Logarithm of malaria incidence rate per 10,000") +
  theme_classic()


mu3 <- as.matrix(cbind(1, d3$t)) %*% 
  as.matrix(fit3$summary.fixed[1:2, ])

nu3 <- as.matrix(cbind(d3$z, d3$s)) %*% 
  as.matrix(fit3$summary.fixed[3:4, ])
nu3[nu3 == 0] <- NA
nu3 <- mu3 + nu3
mu3[!is.na(nu3)] <- NA

g3 <- d3 %>% 
  ggplot(aes(x=date, y=log(y))) + 
  geom_line() +
  geom_line(aes(x=date, y=fit3$summary.fitted.values$mean), col="green4", linetype=2) +
  geom_ribbon(aes(ymin=fit3$summary.fitted.values$`0.025quant`,
                  ymax=fit3$summary.fitted.values$`0.975quant`),
              col="green", fill="green", alpha=0.5) + 
  geom_line(aes(x=date, y=mu3[, 1]), col="blue") +
  geom_ribbon(aes(ymin=mu3[, 3], ymax=mu3[, 5]), col="blue", fill="blue", alpha=0.5) +
  geom_line(aes(x=date, y=nu3[, 1]), col="red") +
  geom_ribbon(aes(ymin=nu3[, 3], ymax=nu3[, 5]), col="red", fill="red", alpha=0.5) +
  labs(x="Date", y="Logarithm of malaria incidence rate per 10,000") +
  theme_classic()

ggarrange(g1, g2, g3, nrow=3, ncol=1,
          labels=c("Batu", "Dire Dawa", "Metehara"),
          label.x=0.05)

library("devEMF")
ggsave(filename="Fig2.emf", device=emf, width=10, height=11)


