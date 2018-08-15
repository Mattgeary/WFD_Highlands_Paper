

## Load libraries
library(unmarked)
library(raster)
library(AICcmodavg)

## Read in occupancy data
WFD.occu.mat <- read.csv("Data/WFDHighland_wide_UMF_FULL.csv")[,2:17]

### Aggregate occupancy matrix to 5 year periods
WFD.occu.5.mat <- data.frame("p1" = apply(WFD.occu.mat[,1:5], 1, max, na.rm = TRUE), "p2" = apply(WFD.occu.mat[,6:10], 1, max, na.rm = TRUE), "p3" = apply(WFD.occu.mat[,11:16], 1, max, na.rm = TRUE))
WFD.occu.5.mat[WFD.occu.5.mat == -Inf] <- NA

a <- which(!is.na(apply(WFD.occu.5.mat[,1:2], 1, max)))
a <- append(a, !is.na(apply(WFD.occu.5.mat[,1:3], 1, max)))
a <- append(a, !is.na(apply(WFD.occu.5.mat[,2:3], 1, max)))
a <- unique(a)

WFD.occu.a <- WFD.occu.5.mat[a,]

# Read in site and observation covariates
WFD.siteCovs <- read.csv("Data/EnvVars_FULL.csv")
WFD.siteCovs.a <- WFD.siteCovs[a,]

## Identify rows with NA values - these are outside of the Scottish mainland
NA.vals <- which(is.na(WFD.siteCovs.a$AWI))
NA.vals <- append(NA.vals, which(is.na(WFD.siteCovs.a$pet_emerge)))
NA.vals <- append(NA.vals, which(is.na(WFD.siteCovs.a$wildness)))
NA.vals <- append(NA.vals, which(is.na(WFD.siteCovs.a$peat)))
NA.vals <- unique(NA.vals)

WFD.obsCov.wind <- read.csv("Data/WFDObsCov_wind.csv")
WFD.obsCov.wind.5 <- data.frame("p1" = apply(WFD.obsCov.wind[,1:5], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.wind[,6:10], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.wind[,11:16], 1, mean, na.rm = TRUE))
WFD.obsCov.wind.5[WFD.obsCov.wind.5 == -Inf] <- NA
WFD.obsCov.wind.a <- WFD.obsCov.wind.5[a,]
WFD.obsCov.tas <- read.csv("Data/WFDObsCov_tas.csv")
WFD.obsCov.tas.5 <- data.frame("p1" = apply(WFD.obsCov.tas[,1:5], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.tas[,6:10], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.tas[,11:16], 1, mean, na.rm = TRUE))
WFD.obsCov.tas.5[WFD.obsCov.tas.5 == -Inf] <- NA
WFD.obsCov.tas.a <- WFD.obsCov.tas.5[a,]

WFD.occu.a <- WFD.occu.a[-c(NA.vals),]
WFD.siteCovs.a <- WFD.siteCovs.a[-c(NA.vals),]
WFD.obsCov.tas.a <- WFD.obsCov.tas.a[-c(NA.vals),]
WFD.obsCov.wind.a <- WFD.obsCov.wind.a[-c(NA.vals),]

# Missing covariates indicated by model fit below
#WFD.occu.mat <- WFD.occu.mat[-c(1016, 1027, 1037),]
#WFD.siteCovs <- WFD.siteCovs[-c(1016, 1027, 1037),]
#WFD.obsCov.tas.5 <- WFD.obsCov.tas.5[-c(1016, 1027, 1037),]
#WFD.obsCov.wind.5 <- WFD.obsCov.wind.5[-c(1016, 1027, 1037),]


WFD.obsCov.AWI <- as.data.frame(WFD.siteCovs.a$AWI)
for(i in 2:3){
WFD.obsCov.AWI <- cbind(WFD.obsCov.AWI, as.data.frame(WFD.siteCovs.a$AWI))
}
WFD.obsCov.wildness <- as.data.frame(WFD.siteCovs.a$wildness)
for(i in 2:3){
WFD.obsCov.wildness <- cbind(WFD.obsCov.wildness, as.data.frame(WFD.siteCovs.a$wildness))
}

## Create observation covariates list
WFD.obsCov <- list('AWI' = WFD.obsCov.AWI, 'wildness' = WFD.obsCov.wildness, 'wind' = WFD.obsCov.wind.a, 'tas' = WFD.obsCov.tas.a)

## Combine conifer/pine and bog/mire/bogwoodland
WFD.siteCovs.a$coniferpine <- WFD.siteCovs.a$conifer + WFD.siteCovs.a$pineforest
WFD.siteCovs.a$wet <- WFD.siteCovs.a$bogs + WFD.siteCovs.a$bogwoodland + WFD.siteCovs.a$mire

# Create Unmarked data frame
WFD.UMF <- unmarkedFrameOccu(y = WFD.occu.a, siteCovs = WFD.siteCovs.a, obsCovs = WFD.obsCov)

# Occupancy model

### Candidate models obsCovs
WFD.occu.obs.1 <- occu(~1 ~1, WFD.UMF)
WFD.occu.obs.2 <- occu(~AWI ~1, WFD.UMF)
WFD.occu.obs.3 <- occu(~wildness ~1, WFD.UMF)
WFD.occu.obs.4 <- occu(~wind ~1, WFD.UMF)
#WFD.occu.obs.5 <- occu(~tas ~1, WFD.UMF)
#WFD.occu.obs.6 <- occu(~AWI + wildness + wind + tas ~1, WFD.UMF)
WFD.occu.obs.7 <- occu(~AWI + wildness + wind ~1, WFD.UMF)
WFD.occu.obs.8 <- occu(~AWI + wildness + tas ~1, WFD.UMF)
WFD.occu.obs.9 <- occu(~AWI + wildness + tas ~1, WFD.UMF)
#WFD.occu.obs.10 <- occu(~AWI + wind + tas ~1, WFD.UMF)
#WFD.occu.obs.11 <- occu(~wildness + wind + tas ~1, WFD.UMF)
WFD.occu.obs.12 <- occu(~AWI + wildness ~1, WFD.UMF)
WFD.occu.obs.13 <- occu(~AWI + wind ~1, WFD.UMF)
WFD.occu.obs.14 <- occu(~wildness + wind ~1, WFD.UMF)
#WFD.occu.obs.15 <- occu(~AWI + tas ~1, WFD.UMF)
WFD.occu.obs.16 <- occu(~wildness + tas ~1, WFD.UMF)
#WFD.occu.obs.17 <- occu(~wind + tas ~1, WFD.UMF)

#WFD.obs.fits <- fitList(WFD.occu.obs.1, WFD.occu.obs.2, WFD.occu.obs.3, WFD.occu.obs.4, WFD.occu.obs.5, WFD.occu.obs.6,  WFD.occu.obs.7, WFD.occu.obs.8, WFD.occu.obs.9, WFD.occu.obs.10, WFD.occu.obs.11, WFD.occu.obs.12, WFD.occu.obs.13, WFD.occu.obs.14, WFD.occu.obs.15, WFD.occu.obs.16, WFD.occu.obs.17)
WFD.obs.fits <- fitList(WFD.occu.obs.1, WFD.occu.obs.2, WFD.occu.obs.3, WFD.occu.obs.4, WFD.occu.obs.7, WFD.occu.obs.8, WFD.occu.obs.9, WFD.occu.obs.12, WFD.occu.obs.13, WFD.occu.obs.14, WFD.occu.obs.16)
WFD.modSel.obs <- modSel(WFD.obs.fits)
WFD.modSel.obs



### Candidate models siteCovs
WFD.occu.1 <- occu(~wildness + wind ~1, WFD.UMF)
WFD.occu.2 <- occu(~wildness + wind~ pineforest*wet + pet_emerge, WFD.UMF)
WFD.occu.3 <- occu(~wildness + wind ~pineforest + pet_emerge + bogs, WFD.UMF)
WFD.occu.4 <- occu(~wildness + wind ~pineforest + pet_emerge + wet, WFD.UMF)
WFD.occu.5 <- occu(~wildness + wind ~pineforest + bogs, WFD.UMF)
WFD.occu.6 <- occu(~wildness + wind ~pineforest + wet, WFD.UMF)
WFD.occu.7 <- occu(~wildness + wind ~pineforest + pet_emerge, WFD.UMF)
WFD.occu.8 <- occu(~wildness + wind ~pineforest*pet_emerge + wet, WFD.UMF)
WFD.occu.9 <- occu(~wildness + wind ~pineforest*pet_emerge + bogs, WFD.UMF)
WFD.occu.10 <- occu(~wildness + wind ~pineforest*pet_emerge, WFD.UMF)
WFD.occu.11 <- occu(~wildness + wind ~pineforest*wet, WFD.UMF)
WFD.occu.12 <- occu(~wildness + wind ~pineforest*bogs, WFD.UMF)
WFD.occu.13 <- occu(~wildness + wind ~pineforest, WFD.UMF)
WFD.occu.14 <- occu(~wildness + wind ~coniferpine + pet_emerge + wet, WFD.UMF)
## Models below don't converge
WFD.occu.15 <- occu(~wildness + wind ~coniferpine*wet + pet_emerge, WFD.UMF)
WFD.occu.16 <- occu(~wildness + wind ~coniferpine*bogs + pet_emerge, WFD.UMF)

WFD.site.fits <- fitList(WFD.occu.1, WFD.occu.2, WFD.occu.3, WFD.occu.4, WFD.occu.5, WFD.occu.6, WFD.occu.7, WFD.occu.8, WFD.occu.9, WFD.occu.10, WFD.occu.11, WFD.occu.12, WFD.occu.13, WFD.occu.14, WFD.occu.15, WFD.occu.16)
#WFD.site.fits <- fitList(WFD.occu.1, WFD.occu.2, WFD.occu.3, WFD.occu.4, WFD.occu.5, WFD.occu.6, WFD.occu.7, WFD.occu.8, WFD.occu.9, WFD.occu.10, WFD.occu.11, WFD.occu.12, WFD.occu.13, WFD.occu.14)
WFD.modSel <- modSel(WFD.site.fits)
WFD.modSel

WFD.gof <- mb.gof.test(WFD.occu.3, nsim = 1000)

