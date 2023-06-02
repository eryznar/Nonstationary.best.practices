# INSTALL PACKAGES --------------------------------------------------------------
#install.packages(c("gtools", "ncdf4", "zoo", "scales", "nlme", "gplots", "maps", "mapdata", "chron", "fields",
#"pracma", "FactoMineR", "lmtest", "MuMIn", "broom", "reshape2"))

# LOAD PACKAGES -----------------------------------------------------------------
library(gtools)
library(ncdf4)
library(zoo)
library(scales) 
library(nlme)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(pracma)
library(FactoMineR)
library(lmtest)
library(MuMIn)
library(broom)
library(reshape2)
library(ggplot2)
library(purrr)

# LOAD DATA ---------------------------------------------------------------------
nc <- nc_open("./data/nceiErsstv5_4ec2_7c0e_3c10.nc")

# PROCESS SST DATA --------------------------------------------------------------
# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
sst.d <- dates(h, origin = c(1,1,1970))

sst.x <- ncvar_get(nc, "longitude")
sst.y <- ncvar_get(nc, "latitude")

# save months and years for use later on
m <- months(sst.d)
yrs <- years(sst.d)

# and set a couple functions for standardizing below
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

SST <- ncvar_get(nc,  "sst")
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))

# Keep track of corresponding latitudes and longitudes of each column:
sst.lat <- rep(sst.y, length(sst.x))
sst.lon <- rep(sst.x, each = length(sst.y))
dimnames(SST) <- list(as.character(sst.d), paste("N", sst.lat, "E", sst.lon, sep=""))

# Extract study area
# 54-62N and 200-226E
keep <- sst.lat %in% 54:62 & sst.lon %in% 200:226
SST <- SST[,keep]
sst.lat <- sst.lat[keep]
sst.lon <- sst.lon[keep]
sst.y <- sst.y[sst.y %in% 54:62]
sst.x <- sst.x[sst.x %in% 200:226]

# drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200")
SST[,BB] <- NA

# now get monthly means weighted by area, using an arithmetic mean
weight <- sqrt(cos(sst.lat*pi/180))
SST.mu <- apply(SST, 1, function(x) weighted.mean(x, weight, na.rm=T))

# now separate out winter
win <- c("Nov", "Dec", "Jan", "Feb", "Mar")
mean <- data.frame(year=yrs, month=m, mean=SST.mu)
mean$win.yr <- as.numeric(as.character(mean$year))
# assigns Nov and Dec to "winter year"" corresponding to Jan
mean$win.yr[mean$month %in% c("Nov", "Dec")] <- mean$win.yr[mean$month %in% c("Nov", "Dec")] + 1  
win.mean <- mean[mean$month %in% win,] # separate out NDJFM
SST.win <- tapply(win.mean$mean, win.mean$win.yr, mean)

# LOAD AND PROCESS CATCH DATA (MIKE'S CODE) --------------------------------------
#Now define the best timing of the change. First step is to calculate PCA of catch.
data <- read.csv("./Data/total_goa_catch.csv", row.names=1)
data$year <- 1900:2014
pc.dat <- data.frame(sock=scale(data$Sockeye[data$year %in% 1965:2012]),
                     pink=scale(data$Pink[data$year %in% 1965:2012]),
                     chum=scale(data$Chum[data$year %in% 1965:2012]),
                     coho=scale(data$Coho[data$year %in% 1965:2012]))
pca <- prcomp(pc.dat)

#These are the loadings and proportion of variance for catch PCA.
pca$rotation
summary(pca)

#Now the model selection process comparing different threshold years. The explanatory variable is GOA winter SST, smoothed with a 3-yr rolling mean.
temp <- SST.win[names(SST.win) %in% 1900:2014]
sst.3 <- rollapply(temp, 3, mean, na.rm=T, fill=NA) 
names(sst.3) <- 1900:2014

dat.thr <- data.frame(pc1=pca$x[,1], sst.3=sst.3[names(sst.3) %in% 1965:2012], year=1965:2012)
thr <- 1974:2002 # these candidate thresholds maintain at least 20% of the data in each era
pc.out <- NA # object to hold results (AICc scores)
for(i in 1:length(thr)){ # loop each threshold
  # i <-1
  dat.thr$era <- "early"
  dat.thr$era[dat.thr$year > thr[i]] <- "late"
  mod <- gls(pc1 ~ sst.3*era, data=dat.thr, correlation = corAR1(form = ~1), method="ML")
  pc.out[i] <- AICc(mod)
}

#Now compare the sst-catch regression pre/post 88/89
dat <- data.frame(pc1=pca$x[,1], sst=sst.3[names(sst.3) %in% 1965:2012], 
                  year=1965:2012)
dat$era <- "early"
dat$era[dat$year > 1988] <- "late"


#And now the era-specific p-values for catch PC1 on SST. Again, these use GLS models allowing for autocorrelated residuals.
dat1 <- dat[dat$era=="early",]
dat2 <- dat[dat$era=="late",]


# SIMULATING TIMESERIES ---------------------------------------------------------

# Specify data frame to hold function outputs
out <- as.data.frame(matrix(NA, nrow = length(ar), ncol = 10))
colnames(out) <- c("iter", "thr", "ar", "ar1", "AIC.1", "AIC.2", "dAIC", "dAIC.thr", "p_stat", "p_nstat")


# Create function to simulate ts and record function outputs
thr.fun <- function(dat.thr, thr){
  dat.thr$era <- "early"
  dat.thr$era[dat.thr$year > thr] <- "late"
  
  return(cbind(thr, dat.thr))
}

thresh %>%
  map_df(~thr.fun(dat.thr, .x)) -> dat.thr

dat.thr <- data.frame(pc1=pca$x[,1], sst.3=sst.3[names(sst.3) %in% 1965:2012], year=1965:2012)
thr <- 1974:2002 # these candidate thresholds maintain at least 20% of the data in each era




ar <- seq(0.2, 0.9, by = 0.1)
iter <- seq(1, 1000, by = 1)

dat <- rbind(dat1, dat2)

out <- as.data.frame(matrix(NA, nrow = length(ar) * length(thr), ncol = 9))
colnames(out) <- c("thr", "ar", "ar1", "AIC.1", "AIC.2", "dAIC", "dAIC.thr", "p_stat", "p_nstat")

expand.grid(thr = thr, ar = ar) -> cc

for(ii in 1:length(ar)){
  sim.fun <- function(dat, dat.thr, thr, ar){
    for(ii in 1:length(thr)){
      
      #Generate random ts
      ts.sim <- arima.sim(model = list(order = c(1,0,0), ar = ar), n = 48)
      
      #Calculate autocorrelation at lag 1
      ar1 <- acf(ts.sim)$acf[2]
      
      #Generate moving threshold data
      dat.thr$era <- "early"
      dat.thr$era[dat.thr$year > thr[ii]] <- "late"
      
      #Bind simulated ts with sst and era data
      dat %>%
        select(!pc1) %>%
        cbind(., ts.sim) -> sim.dat
      
      dat.thr %>%
        select(!pc1) %>%
        cbind(.,ts.sim) -> sim.dat.thr
      
      
      #Run models with both terms
      mod1 <- gls(ts.sim ~ sst, data = sim.dat, correlation = corAR1())
      mod2 <- gls(ts.sim ~ sst*era, data = sim.dat, correlation = corAR1())
      
      mod1.thr <- gls(ts.sim ~ sst.3, data = sim.dat.thr, correlation = corAR1())
      mod2.thr <- gls(ts.sim ~ sst.3*era, data = sim.dat.thr, correlation = corAR1())
      
      #Calculate model AICs and difference between
      AIC.1 <- AICc(mod1)
      AIC.2 <- AICc(mod2)
      
      dAIC <- AIC.1 - AIC.2
      
      AIC.1.thr <- AICc(mod1.thr)
      AIC.2.thr <- AICc(mod2.thr)
      
      dAIC.thr <- AIC.1.thr - AIC.2.thr
      
      p_stat <- summary(mod1)$tTable[2,4]
      p_nstat <- summary(mod2)$tTable[2,4]
      
    }
    out <- data.frame(thr = thr[ii], ar = ar, ar1 = ar1,AIC.1 = AIC.1, AIC.2 = AIC.2, 
                      dAIC = dAIC, dAIC.thr = dAIC.thr, p_stat = p_stat, p_nstat = p_nstat)
    return(out)
  }
  
  
  
  
  ar %>%
    map_df(sim.fun(dat, dat.thr, thr, ar))-> jj
  
  map2_df(cc$thr, cc$ar, sim.fun(dat, dat.thr, cc$thr, cc$ar))
  
  out <- as.data.frame(matrix(NA, nrow = nrow(cc), ncol = 10))
  colnames(out) <- c("iter", "thr", "ar", "ar1", "AIC.1", "AIC.2", "dAIC", "dAIC.thr", "p_stat", "p_nstat")
  
  sim <- function(cc, iter){
    rep(
      for(ii in 1:nrow(cc)){
        
        #Generate random ts
        ts.sim <- arima.sim(model = list(order = c(1,0,0), ar = cc$ar[ii]), n = 48)
        
        #Calculate autocorrelation at lag 1
        ar1 <- acf(ts.sim)$acf[2]
        
        #Generate moving threshold data
        dat.thr$era <- "early"
        dat.thr$era[dat.thr$year > cc$thr[ii]] <- "late"
        
        #Bind simulated ts with sst and era data
        dat %>%
          select(!pc1) %>%
          cbind(., ts.sim) -> sim.dat
        
        dat.thr %>%
          select(!pc1) %>%
          cbind(.,ts.sim) -> sim.dat.thr
        
        
        #Run models with both terms
        mod1 <- gls(ts.sim ~ sst, data = sim.dat, correlation = corAR1())
        mod2 <- gls(ts.sim ~ sst*era, data = sim.dat, correlation = corAR1())
        
        mod1.thr <- gls(ts.sim ~ sst.3, data = sim.dat.thr, correlation = corAR1())
        mod2.thr <- gls(ts.sim ~ sst.3*era, data = sim.dat.thr, correlation = corAR1())
        
        #Calculate model AICs and difference between
        AIC.1 <- AICc(mod1)
        AIC.2 <- AICc(mod2)
        
        dAIC <- AIC.1 - AIC.2
        
        AIC.1.thr <- AICc(mod1.thr)
        AIC.2.thr <- AICc(mod2.thr)
        
        dAIC.thr <- AIC.1.thr - AIC.2.thr
        
        p_stat <- summary(mod1)$tTable[2,4]
        p_nstat <- summary(mod2)$tTable[2,4]
        
        out[ii,] <- data.frame(iter = iter, thr = cc$thr[ii], ar = cc$ar[ii], ar1 = ar1,AIC.1 = AIC.1, AIC.2 = AIC.2, 
                               dAIC = dAIC, dAIC.thr = dAIC.thr, p_stat = p_stat, p_nstat = p_nstat)
      }
      ,iter)
    
    return(out)
  }
  
  
  iter %>%
    map_df(~sim(cc, .x)) -> jj
  
  sim.out <- jj
  # Bin data, calculate proportion of dAIC values >= 1 in each bin
  sim.out %>%
    mutate(bin = case_when((abs(ar1) <= 0.1) ~ "≤0.1",
                           (abs(ar1) > 0.1 & abs(ar1) <= 0.2) ~ "0.11-0.2",
                           (abs(ar1) > 0.2 & abs(ar1) <= 0.3) ~ "0.21-0.3",
                           (abs(ar1) > 0.3 & abs(ar1) <= 0.4) ~ "0.31-0.4",
                           (abs(ar1) > 0.4 & abs(ar1) <= 0.5) ~ "0.41-0.5",
                           (abs(ar1) > 0.5 & abs(ar1) <= 0.6) ~ "0.51-0.6",
                           (abs(ar1) > 0.6 & abs(ar1) <= 0.7) ~ "0.61-0.7",
                           (abs(ar1) > 0.7 & abs(ar1) <= 0.8) ~ "0.71-0.8",
                           (abs(ar1) > 0.8 & abs(ar1) <= 0.9) ~ "0.81-0.9",
                           (abs(ar1) > 0.9 & abs(ar1) <= 1) ~ "0.91-1")) -> binned.dat
  
  binned.dat %>%
    group_by(bin) %>%
    reframe(N = n(),
            n_pos = sum(dAIC >=1),
            n_pos.thr = sum(dAIC.thr >=1),
            prop_pos = n_pos/N,
            prop_pos.thr = n_pos.thr/N) -> plot.dat
  
  binned.dat %>%
    melt(c("p_stat", "p_nstat"), id.vars = "bin") %>%
    mutate(variable = ifelse(variable == "p_stat", "Stationary", "Non-stationary")) %>%
    group_by(bin, variable) %>%
    reframe(N = n(),
            n_sig = sum(value <0.05),
            prop_sig = n_sig/N) -> p.dat
  
  
  
  # Plot
  ggplot(p.dat, aes(x = bin, y = prop_sig, color = variable, group = variable)) +
    geom_line(size = 1.25) +
    geom_point(size = 2.5)+
    scale_color_manual(values = c("#C0C9FE", "#ABD89D"), labels = c("Non-stationary model", "Stationary model"),
                       name = NULL)+
    theme_bw() +
    theme(legend.position = "bottom")+
    xlab("AR(1)")+
    ylab("Proportion p<0.05") -> pval_plot
  
  ggsave(plot = pval_plot, "./Figures/pval_plot.png", 
         height=5, width=7, units="in")
  
  ggplot(plot.dat) +
    geom_line(aes(x = bin, y = prop_pos, color = "#C0C9FE"), group = 1, linewidth = 1.25)+
    geom_line(aes(x = bin, y = prop_pos.thr, , color = "#ABD89D"), group = 1, linewidth = 1.25)+
    geom_point(aes(x = bin, y = prop_pos, color = "#C0C9FE"), size = 2.5) +
    geom_point(aes(x = bin, y = prop_pos.thr, , color = "#ABD89D"), size = 2.5)+
    scale_color_manual(values = c("#C0C9FE", "#ABD89D"), labels = c("Stationary threshold", "Moving threshold"),
                       name = NULL)+
    theme_bw() +
    theme(legend.position =  "bottom")+
    xlab("AR(1)")+
    ylab("dAIC \u2265 1") -> dAIC_plot
  
  ggsave(plot = dAIC_plot, "./Figures/dAIC_plot.png", 
         height=5, width=7, units="in")
  