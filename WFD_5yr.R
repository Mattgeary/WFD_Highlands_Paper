#### Occupancy model using summarised occupancy matrix 5 yr, 5 yr, 6 yr (bit of a fudge). Seems to work better in some ways but the 'best' models don't converge

## Load libraries
library(unmarked)
library(raster)
library(AICcmodavg)

## Read in occupancy data
WFD.occu.mat <- read.csv("Data/WFDHighland_wide_UMF_FULL.csv")[,2:17]

### Aggregate occupancy matrix to 5 year periods
WFD.occu.5.mat <- data.frame("p1" = apply(WFD.occu.mat[,1:5], 1, max, na.rm = TRUE), "p3" = apply(WFD.occu.mat[,6:10], 1, max, na.rm = TRUE), "p3" = apply(WFD.occu.mat[,11:16], 1, max, na.rm = TRUE))
WFD.occu.5.mat[WFD.occu.5.mat == -Inf] <- NA


# Read in site and observation covariates
WFD.siteCovs <- read.csv("Data/EnvVars_FULL.csv")


## Identify rows with NA values - these are outside of the Scottish mainland
NA.vals <- which(is.na(WFD.siteCovs$AWI))
NA.vals <- append(NA.vals, which(is.na(WFD.siteCovs$pet_emerge)))
NA.vals <- append(NA.vals, which(is.na(WFD.siteCovs$wildness)))
NA.vals <- append(NA.vals, which(is.na(WFD.siteCovs$peat)))

WFD.obsCov.wind <- read.csv("Data/WFDObsCov_wind.csv")
WFD.obsCov.wind.5 <- data.frame("p1" = apply(WFD.obsCov.wind[,1:5], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.wind[,6:10], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.wind[,11:16], 1, mean, na.rm = TRUE))
WFD.obsCov.wind.5[WFD.obsCov.wind.5 == -Inf] <- NA
WFD.obsCov.tas <- read.csv("Data/WFDObsCov_tas.csv")
WFD.obsCov.tas.5 <- data.frame("p1" = apply(WFD.obsCov.tas[,1:5], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.tas[,6:10], 1, mean, na.rm = TRUE), "p3" = apply(WFD.obsCov.tas[,11:16], 1, mean, na.rm = TRUE))
WFD.obsCov.tas.5[WFD.obsCov.tas.5 == -Inf] <- NA


WFD.occu.mat <- WFD.occu.5.mat[-c(NA.vals),]
WFD.siteCovs <- WFD.siteCovs[-c(NA.vals),]
WFD.obsCov.tas.5 <- WFD.obsCov.tas.5[-c(NA.vals),]
WFD.obsCov.wind.5 <- WFD.obsCov.wind.5[-c(NA.vals),]

# Missing covariates indicated by model fit below
WFD.occu.mat <- WFD.occu.mat[-c(1016, 1027, 1037),]
WFD.siteCovs <- WFD.siteCovs[-c(1016, 1027, 1037),]
WFD.obsCov.tas.5 <- WFD.obsCov.tas.5[-c(1016, 1027, 1037),]
WFD.obsCov.wind.5 <- WFD.obsCov.wind.5[-c(1016, 1027, 1037),]


WFD.obsCov.AWI <- as.data.frame(WFD.siteCovs$AWI)
for(i in 2:3){
WFD.obsCov.AWI <- cbind(WFD.obsCov.AWI, as.data.frame(WFD.siteCovs$AWI))
}
WFD.obsCov.wildness <- as.data.frame(WFD.siteCovs$wildness)
for(i in 2:3){
WFD.obsCov.wildness <- cbind(WFD.obsCov.wildness, as.data.frame(WFD.siteCovs$wildness))
}

## Create observation covariates list
WFD.obsCov <- list('AWI' = WFD.obsCov.AWI, 'wildness' = WFD.obsCov.wildness, 'wind' = WFD.obsCov.wind.5, 'tas' = WFD.obsCov.tas.5)

## Combine conifer/pine and bog/mire/bogwoodland
WFD.siteCovs$coniferpine <- WFD.siteCovs$conifer + WFD.siteCovs$pineforest
WFD.siteCovs$wet <- WFD.siteCovs$bogs + WFD.siteCovs$bogwoodland + WFD.siteCovs$mire



# Create Unmarked data frame
WFD.UMF <- unmarkedFrameOccu(y = WFD.occu.mat, siteCovs = WFD.siteCovs, obsCovs = WFD.obsCov)

# Scale and centre siteCovs
siteCovs(WFD.UMF) <- scale(siteCovs(WFD.UMF))
sC.scale <- scale(siteCovs(WFD.UMF))
obsCovs(WFD.UMF) <- scale(obsCovs(WFD.UMF))
oC.scale <- scale(obsCovs(WFD.UMF))

# Occupancy model

### Candidate models obsCovs
WFD.occu.obs.1 <- occu(~1 ~1, WFD.UMF)
WFD.occu.obs.2 <- occu(~AWI ~1, WFD.UMF)
WFD.occu.obs.3 <- occu(~wildness ~1, WFD.UMF)
WFD.occu.obs.4 <- occu(~wind ~1, WFD.UMF)
WFD.occu.obs.5 <- occu(~tas ~1, WFD.UMF)
WFD.occu.obs.6 <- occu(~AWI + wildness + wind + tas ~1, WFD.UMF)
WFD.occu.obs.7 <- occu(~AWI + wildness + wind ~1, WFD.UMF)
WFD.occu.obs.8 <- occu(~AWI + wildness + tas ~1, WFD.UMF)
WFD.occu.obs.9 <- occu(~AWI + wildness + tas ~1, WFD.UMF)
WFD.occu.obs.10 <- occu(~AWI + wind + tas ~1, WFD.UMF)
WFD.occu.obs.11 <- occu(~wildness + wind + tas ~1, WFD.UMF)
WFD.occu.obs.12 <- occu(~AWI + wildness ~1, WFD.UMF)
WFD.occu.obs.13 <- occu(~AWI + wind ~1, WFD.UMF)
WFD.occu.obs.14 <- occu(~wildness + wind ~1, WFD.UMF)
WFD.occu.obs.15 <- occu(~AWI + tas ~1, WFD.UMF)
WFD.occu.obs.16 <- occu(~wildness + tas ~1, WFD.UMF)
WFD.occu.obs.17 <- occu(~wind + tas ~1, WFD.UMF)

WFD.obs.fits <- fitList(WFD.occu.obs.1, WFD.occu.obs.2, WFD.occu.obs.3, WFD.occu.obs.4, WFD.occu.obs.5, WFD.occu.obs.6,  WFD.occu.obs.7, WFD.occu.obs.8, WFD.occu.obs.9, WFD.occu.obs.10, WFD.occu.obs.11, WFD.occu.obs.12, WFD.occu.obs.13, WFD.occu.obs.14, WFD.occu.obs.15, WFD.occu.obs.16, WFD.occu.obs.17)
WFD.modSel.obs <- modSel(WFD.obs.fits)
WFD.modSel.obs



### Candidate models siteCovs
WFD.occu.1 <- occu(~wildness + wind + tas ~1, WFD.UMF)
WFD.occu.2 <- occu(~wildness + wind + tas~ pineforest*wet + pet_emerge, WFD.UMF)
WFD.occu.3 <- occu(~wildness + wind + tas ~pineforest + pet_emerge + bogs, WFD.UMF)
WFD.occu.4 <- occu(~wildness + wind + tas ~pineforest + pet_emerge + wet, WFD.UMF)
WFD.occu.5 <- occu(~wildness + wind + tas ~pineforest + bogs, WFD.UMF)
WFD.occu.6 <- occu(~wildness + wind + tas ~pineforest + wet, WFD.UMF)
WFD.occu.7 <- occu(~wildness + wind + tas ~pineforest + pet_emerge, WFD.UMF)
WFD.occu.8 <- occu(~wildness + wind + tas ~pineforest*pet_emerge + wet, WFD.UMF)
WFD.occu.9 <- occu(~wildness + wind + tas ~pineforest*pet_emerge + bogs, WFD.UMF)
WFD.occu.10 <- occu(~wildness + wind + tas ~pineforest*pet_emerge, WFD.UMF)
WFD.occu.11 <- occu(~wildness + wind + tas ~pineforest*wet, WFD.UMF)
WFD.occu.12 <- occu(~wildness + wind + tas ~pineforest*bogs, WFD.UMF)
WFD.occu.13 <- occu(~wildness + wind + tas ~pineforest, WFD.UMF)
WFD.occu.14 <- occu(~wildness + wind + tas ~coniferpine + pet_emerge + wet, WFD.UMF)
WFD.occu.15 <- occu(~wildness + wind + tas ~coniferpine*wet + pet_emerge, WFD.UMF)
WFD.occu.16 <- occu(~wildness + wind + tas ~coniferpine*bogs + pet_emerge, WFD.UMF)

WFD.site.fits <- fitList(WFD.occu.1, WFD.occu.2, WFD.occu.3, WFD.occu.4, WFD.occu.5, WFD.occu.6, WFD.occu.7, WFD.occu.8, WFD.occu.9, WFD.occu.10, WFD.occu.11, WFD.occu.12, WFD.occu.13, WFD.occu.14, WFD.occu.15, WFD.occu.16)
WFD.modSel <- modSel(WFD.site.fits)
WFD.modSel
