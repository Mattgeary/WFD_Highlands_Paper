##Example code to perform join count chi squared test on simulated occupancy dataset##
##covariates are allowed in functions and different models can be fit##
##Markovian models which allow detection probability to depend on whether the previous##
##visit detected the species or not can also be fit##
##Also performs the goodness of fit test on the simulated dataset##

##Author: Wilson Wright##

##set working directory to where the functions.R script is saved##

setwd("YourDirectoryHere")

##Load required packages##
##Install all of these packages prior to running code if they are not already installed##
##Also need to install Rtools (for compiling code to C++ and using Rcpp package)##
require(dplyr)
require(tidyr)
require(compiler)
require(inline)
require(Rcpp)
require(RcppArmadillo)
bdiagMat <- metaSEM::bdiagMat #need to install metaSEM package, but dont require() it

##After the above packages are loaded, the following line will make functions available to use##
##this step will take some time to create all the necessary functions##
source("Data1B.R")

##As an example, simulate a dataset with Markovian correlation to use the join count test on##
n <- 40 #number of sample units
J <- 4 #number of 'groups' or detectors per sample unit
K <- 4 #number of revisits per group

#create two types of sites with different occupancy probabilities in each
#this will be in the format needed for the functions to fit models and perform the tests
#set seed to reproduce this example exactly
set.seed(1213)
site.covariates <- data.frame(site = c(1:n), type = factor(sample(c('A', 'B'), n, replace = TRUE)))

#Specify vector of occupancy probabilities based on 'type'
occ.p <- ifelse(site.covariates$type == 'A', 0.4, 0.8)

#Now simulate true occupancy for these sites
true.z <- rbinom(n, 1, occ.p)

##Create continuous covariate that differs by detector within a sample unit
##This will impact detection probabilities
group.IDs <- rep(1:J, n)
group.covariates <- data.frame(site = rep(1:n, each=J), group = group.IDs, det.predictor = rnorm(length(group.IDs)))

##Create another continuous covariate that differs by visit within each detector
##This will not impact detection probabilites in this simulation
visit.covariates <- list(useless.predictor = data.frame(
  matrix(rnorm(n*J*K), ncol=K, nrow=n*J)
))

##Simulate detection history matrix based on det. probs. specified by intercept and slope with group.covariate
##also simulate a Markovian process where the probability of detection depends on the previous detection or not
detection.history <- matrix(NA, ncol=K, nrow=n*J)
for(i in 1:nrow(detection.history)){
  p0 <- plogis(-1 + 0.5*group.covariates$det.predictor[i])
  p1 <- plogis(-1 + 0.5*group.covariates$det.predictor[i] + 2)
  detection.history[i, 1] <- rbinom(1, 1, true.z[group.covariates$site[i]]*p0/(p0 + (1-p1)))
  for(k in 2:K){
    detection.history[i, k] <- rbinom(1, 1, true.z[group.covariates$site[i]]*(p1*detection.history[i, (k-1)] +
                                                                                p0*(1-detection.history[i, (k-1)])))
  }
}

##Fit a basic model and display the estimates##
##First argument is form of occupancy model, second is form of detection model
fit1 <- fit.occ.model(~type, ~det.predictor, detection.history, site.covariates, group.covariates,
                      visit.covariates, mark.corr=FALSE)
occ.model.summary(fit1)

##Do JC chi2 test, will take some time to do all the permutations
##Need to specify same model as fit above. Also include the model fit, number of simulated datasets, 
##True dataset, and the covariates specified above.
jc.test1 <- jc.chi.test(~type, ~det.predictor, fit1, nsims=500, detection.history, site.covariates, group.covariates,
                        visit.covariates)
jc.test1$p.val1 #p-value from using neighbor definition with all revisits at a detector
jc.test1$p.val2 #p-value from using only adjacent revisits as neighbors
##Both of these p-values are fairly small (0.07 and 0.05),##
##suggesting some evidence that there are more joins that expected under the fit of this basic model.

##Next perform the goodness-of-fit test with this model fit
gof1 <- gof.test(~type, ~det.predictor, fit1, nsims=500, detection.history, site.covariates, group.covariates,
                 visit.covariates)
gof1$p.val #large p-value would indicate no evidence for lack of fit of the model

##Next a model can be fit with a Markov process and then the jc test can be performed again##

##Fit a Markov model and display the estimates##
fit2 <- fit.occ.model(~type, ~det.predictor, detection.history, site.covariates, group.covariates,
                      visit.covariates, mark.corr=TRUE)
occ.model.summary(fit2)

##Do JC chi2 test, will take some time to do all the permutations
jc.test2 <- jc.chi.test(~type, ~det.predictor, fit2, nsims=500, detection.history, site.covariates, group.covariates,
                        visit.covariates)
jc.test2$p.val1 #p-value from using neighbor definition with all revisits at a detector
jc.test2$p.val2 #p-value from using only adjacent revisits as neighbors
##These tests no longer provide any evidence that the number of joins differs from what would be expected##
##p-values are 0.708 and 0.798 respectively##
##This indicates that the Markov model appears adequate for these data##

##The goodness-of-fit test can also be performed using this model fit
gof2 <- gof.test(~type, ~det.predictor, fit2, nsims=500, detection.history, site.covariates, group.covariates,
                 visit.covariates)
gof2$p.val #large p-value would indicate no evidence for lack of fit of the model
