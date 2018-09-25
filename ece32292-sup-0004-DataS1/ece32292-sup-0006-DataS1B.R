##R Script for fitting occupancy models with two forms of replication##
##Allows correlation between consecutive detections from one form of replication##
##User specified covariates and model structure is allowed##
##Estimates are obtained using maximum likelihood via the optim function##

##Need to load the R packages that are listed in the Data1A.R script before running this code##

##Author: Wilson Wright##

##Function to fit a model to occupancy data. Specify a formula for occupancy and detection
##det.matrix is a detection matrix, one row for each grouping unit, columns for serial visits at a group
##NAs are allowed in the detection history matrix, but after NA all subsequent visits must be NA for that group
##Specify a matrix of site.covs with rows=number of sites, columns=number of covariates
##Specify a matrix of group.covs with rows=total number of groups, columns=number of group covariates
##Note that the rows must match those of det.matrix
##Observations and covariates of groups from the same site MUST be in consecutive rows
##Within a site, rows should be ordered by ascending number of serial replicates
##visit.covs is a list of each matrices with dimensions of det.matrix,
##each matrix in the list is a separate covariate. in the list, names given should be covariate names
fit.occ.model <- function(occ.form, det.form, det.matrix, site.covs, group.covs, visit.covs,
                          mark.corr, use.det.vec=FALSE, det.vec=NULL){
  occ.model.matrix <- model.matrix(occ.form, data=site.covs)
  
  det.data.frame <- data.frame(dets=as.vector(t(det.matrix)),
                               site=rep(group.covs$site, each=ncol(det.matrix)),
                               group=rep(group.covs$group, each=ncol(det.matrix)),
                               lapply(visit.covs, function(x){as.vector(t(x))}),
                               first.visit=rep(c(1, rep(0, ncol(det.matrix)-1)), nrow(det.matrix)))
  
  dets.covs.frame <- left_join(det.data.frame, left_join(site.covs, group.covs, by='site'), by=c('site', 'group'))
  dets.covs.frame <- filter(dets.covs.frame, !is.na(dets))
  dets.covs.frame$dets <- if(use.det.vec==FALSE){dets.covs.frame$dets} else{det.vec}
  dets.covs.frame <- mutate(dets.covs.frame, prev.visit=lag(dets, 1, default=0))
  
  det.form2 <- if(mark.corr==TRUE){update.formula(det.form, ~. + prev.visit)} else{det.form}
  
  det.model.matrix <- model.matrix(det.form, dets.covs.frame)
  det.model.matrix2 <- model.matrix(det.form2, dets.covs.frame)
  
  occ.parms <- dimnames(occ.model.matrix)[[2]]
  det.parms <- dimnames(det.model.matrix2)[[2]]
  
  logL.fun <- function(parameters, dets=dets.covs.frame$dets, prev.visit=dets.covs.frame$prev.visit,
                       site=dets.covs.frame$site, first.visit=dets.covs.frame$first.visit,
                       occ.model.matrix=occ.model.matrix, det.model.matrix=det.model.matrix,
                       det.model.matrix2=det.model.matrix2, n.occ.par=length(occ.parms),
                       n.det.par=length(det.parms),
                       mark.corr=mark.corr, occ.parms=occ.parms, det.parms=det.parms){
    betas <- parameters[1:n.occ.par]
    alphas <- parameters[(n.occ.par+1): (n.occ.par+n.det.par)]
    
    n.dets <- tapply(dets, site, sum)
    psi.vec <- plogis(occ.model.matrix%*%betas)
    psi.vec <- tapply(psi.vec, unique(site), mean)
    
    p1s <- if(mark.corr==TRUE){plogis(det.model.matrix%*%alphas[-n.det.par] + alphas[n.det.par])} else{
      plogis(det.model.matrix%*%alphas)
    }
    p0s <- if(mark.corr==TRUE){plogis(det.model.matrix%*%alphas[-n.det.par])} else{p1s}
    
    p.means <- p0s / (p0s + (1-p1s))
    
    probs <- first.visit*(dets*p.means + (1-dets)*(1-p.means)) +
      (1-first.visit)*(dets*(prev.visit*p1s + (1-prev.visit)*p0s) + 
                         (1-dets)*(prev.visit*(1-p1s) + (1-prev.visit)*(1-p0s)))
    
    p.prod.vec <- tapply(probs, site, prod)
    logL <- log(psi.vec*p.prod.vec + (1-psi.vec)*(n.dets==0))
    sum(-1*logL)
  }
  
  output <- optim(c(rep(0, length(c(occ.parms, det.parms)))), logL.fun, method="BFGS", hessian=TRUE,
                  dets=dets.covs.frame$dets, prev.visit=dets.covs.frame$prev.visit,
                  site=dets.covs.frame$site, first.visit=dets.covs.frame$first.visit,
                  occ.model.matrix=occ.model.matrix, det.model.matrix=det.model.matrix,
                  det.model.matrix2=det.model.matrix2, n.occ.par=length(occ.parms), n.det.par=length(det.parms),
                  mark.corr=mark.corr, occ.parms=occ.parms, det.parms=det.parms)
  
  estimates <- output$par
  st.errors <- sqrt(diag(solve(output$hessian)))
  vcov.mat <- solve(output$hessian)
  
  return(list(occ.parms=occ.parms, det.parms=det.parms, estimates=estimates, st.errors=st.errors,
              vcov=vcov.mat, occ.formula=occ.form, det.formula=det.form2,
              occ.model.matrix=occ.model.matrix, det.model.matrix.full=det.model.matrix,
              mark.corr=mark.corr, dets.covs.frame=dets.covs.frame))
}

##Nicely print estimates and standard errors from a model fit with the above function
occ.model.summary <- function(model.fit){
  cat("Model Summary", "\n")
  cat("\n")
  cat("Occupancy Formula =", deparse(model.fit$occ.formula), "\n")
  cat("Coefficient Estimates:", "\n")
  occ.coef <- model.fit$occ.parms
  occ.table <- data.frame(model.fit$estimates[1:length(occ.coef)], model.fit$st.errors[1:length(occ.coef)])
  dimnames(occ.table) <- list(occ.coef, c("Estimate", "St. Error"))
  print(occ.table)
  cat("\n")
  cat("Detection Formula =", deparse(model.fit$det.formula), "\n")
  cat("Coefficient Estimates:", "\n")
  det.coef <- model.fit$det.parms
  det.table <- data.frame(model.fit$estimates[(length(occ.coef)+1):(length(occ.coef)+length(det.coef))],
                          model.fit$st.errors[(length(occ.coef)+1):(length(occ.coef)+length(det.coef))])
  dimnames(det.table) <- list(det.coef, c("Estimate", "St. Error"))
  print(det.table)
}

##Make a function that simulates data based on the estimates from a model fit
simulate.det.hist <- function(model.fit, nsims){
  betas <- model.fit$estimates[1:length(model.fit$occ.parms)]
  alphas <- model.fit$estimates[(length(model.fit$occ.parms)+1):length(model.fit$estimates)]
  
  psis <- plogis(model.fit$occ.model.matrix%*%betas)
  true.z <- rbinom(length(psis), 1, psis)
  z.tab <- tapply(true.z, unique(model.fit$dets.covs.frame$site), mean)
  
  p1s <- if(model.fit$mark.corr==TRUE){plogis(model.fit$det.model.matrix%*%alphas[-length(alphas)] + 
                                                alphas[length(alphas)])} else{
                                                  plogis(model.fit$det.model.matrix%*%alphas)
                                                }
  p0s <- if(model.fit$mark.corr==TRUE){plogis(model.fit$det.model.matrix%*%alphas[-length(alphas)])} else{p1s}
  
  p.means <- p0s / (p0s + (1-p1s))
  
  temp <- model.fit$dets.covs.frame
  
  det.mat <- matrix(NA, nrow=nsims, ncol=nrow(temp))
  for(i in 1:nsims){
    for(j in 1:nrow(temp)){
      det.mat[i,j] <- if(temp$first.visit[j]==1){rbinom(1,1, z.tab[paste(temp$site[j])]*p.means[j])} else{
        rbinom(1, 1, z.tab[paste(temp$site[j])]*(p1s[j]*det.mat[i, (j-1)] + p0s[j]*(1-det.mat[i, (j-1)])))
      }
    }
  }
  return(det.mat)
}

##First create a function to be used later#
##Calculates detection histor probabilities##
calc.hist.prob <- function(x, data, psi.table){
  temp.dets <- rep(as.numeric(x), length(unique(data$site)))
  prev <- rep(lag(as.numeric(x), default=0), length(unique(data$site)))
  probs <- data$first.visit*(temp.dets*data$p.means + (1-temp.dets)*(1-data$p.means)) +
    (1-data$first.visit)*(temp.dets*(prev*data$p1s + 
                                       (1-prev)*data$p0s) + 
                            (1-temp.dets)*(prev*(1-data$p1s) + 
                                             (1-prev)*(1-data$p0s)))
  probs.prod.vec <- tapply(probs, data$site, prod)
  out.probs <- psi.table*probs.prod.vec + (1-psi.table)*(sum(x)==0)
  return(sum(out.probs))
}

##Calculates JC test statistic for both neighbor definitions used in the paper
##Requires model fit using the function above
calc.jc.stat2 <- function(model.fit){
  betas <- model.fit$estimates[1:length(model.fit$occ.parms)]
  alphas <- model.fit$estimates[(length(model.fit$occ.parms)+1):length(model.fit$estimates)]
  
  psis <- plogis(model.fit$occ.model.matrix%*%betas)
  psi.tab <- tapply(psis, unique(model.fit$dets.covs.frame$site), mean)
  
  p1s <- if(model.fit$mark.corr==TRUE){plogis(model.fit$det.model.matrix%*%alphas[-length(alphas)] + 
                                                alphas[length(alphas)])} else{
                                                  plogis(model.fit$det.model.matrix%*%alphas)
                                                }
  p0s <- if(model.fit$mark.corr==TRUE){plogis(model.fit$det.model.matrix%*%alphas[-length(alphas)])} else{p1s}
  p.means <- p0s / (p0s + (1-p1s))
  
  temp.frame <- data.frame(model.fit$dets.covs.frame, p1s=p1s, p0s=p0s, p.means=p.means)
  
  num.grps.frame <- temp.frame %>% group_by(site) %>% summarise(num.grps=length(unique(group)))
  dif.grps <- unique(num.grps.frame$num.grps)
  num.dif.grps <- length(dif.grps)
  
  grp.stats1 <- grp.stats2 <- rep(NA, num.dif.grps)
  for(i in 1:num.dif.grps){
    grp.frame <- filter(temp.frame, site %in%
                          num.grps.frame$site[which(num.grps.frame$num.grps==dif.grps[i])])
    num.visits.frame <- grp.frame %>% group_by(site, group) %>% summarise(num.visits=length(dets))
    spread.visits.frame <- num.visits.frame %>% spread(group, num.visits)
    dif.visits <- unique(spread.visits.frame[ , -1], MARGIN=1)
    num.dif.visits <- nrow(dif.visits)
    
    visit.stats1 <- visit.stats2 <- rep(NA, num.dif.visits)
    for(j in 1:num.dif.visits){
      site.keep <- num.visits.frame %>% group_by(site) %>% 
        summarise(test=identical(num.visits, as.integer(dif.visits[j, ])))
      visit.frame <- filter(grp.frame, site %in% site.keep$site[which(site.keep$test==TRUE)])
      psi.tab.new <- psi.tab[factor(site.keep$site[which(site.keep$test==TRUE)])]
      
      obs.det.mat <- data.frame(select(visit.frame, site, dets), 
                                visit=rep(1:sum(dif.visits[j,]), length(unique(visit.frame$site)))) %>% 
        spread(visit, dets) %>% select(-site) %>% as.matrix()
      t.obs.det.mat <- t(obs.det.mat)
      
      det.hists <- as.matrix(expand.grid(rep(list(c(0,1)), ncol(obs.det.mat))))
      exp.det.hists <- rcpp.exp(det.hists, matrix(visit.frame$p0s, ncol=ncol(obs.det.mat), byrow=T),
                                matrix(visit.frame$p1s, ncol=ncol(obs.det.mat), byrow=T),
                                matrix(visit.frame$p.means, ncol=ncol(obs.det.mat), byrow=T),
                                as.numeric(psi.tab.new),
                                visit.frame$first.visit[1:ncol(obs.det.mat)])
      exp.det.hists[1] <- exp.det.hists[1] + sum(1-as.numeric(psi.tab.new))
      
      t.det.hists <- t(det.hists)
      
      mat1s <- mat2s <- list()
      for(k in 1:length(dif.visits[j ,])){
        mat1s[[k]] <- matrix(1, ncol=as.numeric(dif.visits[j,k]), nrow=as.numeric(dif.visits[j,k]))
        diag(mat1s[[k]]) <- 0
        mat2s[[k]] <- matrix(0, ncol=as.numeric(dif.visits[j,k]), nrow=as.numeric(dif.visits[j,k]))
        delta <- row(mat2s[[k]]) - col(mat2s[[k]])
        mat2s[[k]][abs(delta)==1] <- 1
      }
      w.mat1 <- bdiagMat(mat1s)
      w.mat2 <- bdiagMat(mat2s)
      
      det.hists.BB1 <- rcpp.fact(rcpp.BB(w.mat1, det.hists, t.det.hists))
      det.hists.BB2 <- rcpp.fact(rcpp.BB(w.mat2, det.hists, t.det.hists))
      exp.counts1 <- tapply(exp.det.hists, det.hists.BB1, sum)
      exp.counts2 <- tapply(exp.det.hists, det.hists.BB2, sum)
      
      obs.BB1 <- rcpp.fact2(rcpp.BB(w.mat1, obs.det.mat, t.obs.det.mat), as.numeric(levels(det.hists.BB1)))
      obs.BB2 <- rcpp.fact2(rcpp.BB(w.mat2, obs.det.mat, t.obs.det.mat), as.numeric(levels(det.hists.BB2)))
      
      obs.counts1 <- table(obs.BB1)
      obs.counts2 <- table(obs.BB2)
      
      visit.stats1[j] <- sum( (obs.counts1-exp.counts1)^2 / exp.counts1 )
      visit.stats2[j] <- sum( (obs.counts2-exp.counts2)^2 / exp.counts2 )
    }
    grp.stats1[i] <- sum(visit.stats1)
    grp.stats2[i] <- sum(visit.stats2)
  }
  return(list(test.stat1=sum(grp.stats1), test.stat2=sum(grp.stats2)))
}

##The following creates function in C++ to speed computations##
##These are used by the function above##
rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'

code <- '
mat wmat = as<mat>(wmatin);
mat dets = as<mat>(detsin);
mat tdets = as<mat>(tdetsin);
int m = dets.n_rows;
NumericVector out(m);
for(int i=0; i < m; i++){
out[i]= as_scalar(0.5*(dets.row(i)*wmat*tdets.col(i)));
}
return(wrap(out));
'
rcpp.BB <- cxxfunction(signature(wmatin="numeric", detsin="numeric", tdetsin="numeric"),
                       code, plugin="RcppArmadillo", rcpp_inc)

code3 <- '
NumericVector x(xin);
NumericVector levs = sort_unique(x);
IntegerVector out = match(x, levs);
out.attr("levels") = as<CharacterVector>(levs);
out.attr("class") = "factor";
return(out);
'

rcpp.fact <- cxxfunction(signature(xin="numeric"),
                         code3, plugin="RcppArmadillo", rcpp_inc)

code4 <- '
NumericVector x(xin);
NumericVector levs(levsin);
IntegerVector out = match(x, levs);
out.attr("levels") = as<CharacterVector>(levs);
out.attr("class") = "factor";
return(out);
'

rcpp.fact2 <- cxxfunction(signature(xin="numeric", levsin="numeric"),
                          code4, plugin="RcppArmadillo", rcpp_inc)

code2.new <- '
mat dets = as<mat>(detsin);
mat p0 = as<mat>(p0vec);
mat p1 = as<mat>(p1vec);
mat pm = as<mat>(pmvec);
NumericVector psi(psivec);
NumericVector visit1(firstvisit);

int m = dets.n_rows;
int n = p0.n_rows;
int o = p0.n_cols;
NumericVector out(m);

for(int i=0; i < m; i++){
double x1 = as_scalar(dets(i,0));
NumericVector prev(o);
NumericVector prob(n);
prev[0] = 0;
out[i] = 0;

for(int j=0; j <n; j++){
double temp;
temp = psi[j]*(x1*as_scalar(pm(j, 0)) + (1-x1)*(1-as_scalar(pm(j, 0))));

for(int k=1; k<o; k++){
prev[k] = as_scalar(dets(i, (k-1)));
double x2 = as_scalar(dets(i, k));
temp *= visit1[k]*(x2*as_scalar(pm(j, k)) + (1-x2)*(1-as_scalar(pm(j,k)))) + 
(1-visit1[k])*(x2*(prev[k]*as_scalar(p1(j, k)) + (1-prev[k])*as_scalar(p0(j,k))) + 
(1-x2)*(prev[k]*(1-as_scalar(p1(j, k))) + (1-prev[k])*(1-as_scalar(p0(j,k)))));
}

prob[j] = temp;
out[i] += prob[j];
}
}
return(wrap(out));
'

rcpp.exp <- cxxfunction(signature(detsin="numeric", p0vec="numeric", p1vec="numeric", pmvec = "numeric",
                                  psivec="numeric", firstvisit="numeric"),
                        code2.new, plugin="RcppArmadillo", rcpp_inc)

##The next function performs the test using a specified number of simulated datasets
##Need to specify the model fit (including covariates) and the model object created earlier
jc.chi.test <-function(occ.form, det.form, model.fit, nsims, det.matrix, site.covs, group.covs, visit.covs){
  sim.dets <- simulate.det.hist(model.fit, nsims)
  jc.sims <- apply(sim.dets, 1, function(x){
    calc.jc.stat2(fit.occ.model(occ.form, det.form, det.matrix, site.covs, 
                               group.covs, visit.covs,
                               mark.corr=model.fit$mark.corr, use.det.vec=TRUE, det.vec=x))
  })
  obs.stats <- calc.jc.stat2(model.fit)
  pval1 <- mean(obs.stats[[1]] < unlist(jc.sims)[which(names(unlist(jc.sims))=="test.stat1")]) 
  pval2 <- mean(obs.stats[[2]] < unlist(jc.sims)[which(names(unlist(jc.sims))=="test.stat2")])
  
  return(list(obs.stat1 = obs.stats[[1]], obs.stat2 = obs.stats[[2]],
              p.val1 = pval1, p.val2 = pval2))
}

##GOF test statistic function
calc.gof.stat <- function(model.fit){
  betas <- model.fit$estimates[1:length(model.fit$occ.parms)]
  alphas <- model.fit$estimates[(length(model.fit$occ.parms)+1):length(model.fit$estimates)]
  
  psis <- plogis(model.fit$occ.model.matrix%*%betas)
  psi.tab <- tapply(psis, unique(model.fit$dets.covs.frame$site), mean)
  
  p1s <- if(model.fit$mark.corr==TRUE){plogis(model.fit$det.model.matrix%*%alphas[-length(alphas)] + 
                                                alphas[length(alphas)])} else{
                                                  plogis(model.fit$det.model.matrix%*%alphas)
                                                }
  p0s <- if(model.fit$mark.corr==TRUE){plogis(model.fit$det.model.matrix%*%alphas[-length(alphas)])} else{p1s}
  p.means <- p0s / (p0s + (1-p1s))
  
  temp.frame <- data.frame(model.fit$dets.covs.frame, p1s=p1s, p0s=p0s, p.means=p.means)
  
  num.grps.frame <- temp.frame %>% group_by(site) %>% summarise(num.grps=length(unique(group)))
  dif.grps <- unique(num.grps.frame$num.grps)
  num.dif.grps <- length(dif.grps)
  
  grp.stats <- rep(NA, num.dif.grps)
  for(i in 1:num.dif.grps){
    grp.frame <- filter(temp.frame, site %in%
                          num.grps.frame$site[which(num.grps.frame$num.grps==dif.grps[i])])
    num.visits.frame <- grp.frame %>% group_by(site, group) %>% summarise(num.visits=length(dets))
    spread.visits.frame <- num.visits.frame %>% spread(group, num.visits)
    dif.visits <- unique(spread.visits.frame[ , -1], MARGIN=1)
    num.dif.visits <- nrow(dif.visits)
    
    visit.stats <- rep(NA, num.dif.visits)
    for(j in 1:num.dif.visits){
      site.keep <- num.visits.frame %>% group_by(site) %>% 
        summarise(test=identical(num.visits, as.integer(dif.visits[j, ])))
      visit.frame <- filter(grp.frame, site %in% site.keep$site[which(site.keep$test==TRUE)])
      psi.tab.new <- psi.tab[factor(site.keep$site[which(site.keep$test==TRUE)])]
      
      obs.det.mat <- data.frame(select(visit.frame, site, dets), 
                                visit=rep(1:sum(dif.visits[j,]), length(unique(visit.frame$site)))) %>% 
        spread(visit, dets) %>% select(-site) 
      obs.det.hists <- obs.det.mat %>% unique(MARGIN=1)
      
      obs.hists.counts <- apply(obs.det.hists, 1, FUN=function(x){sum(apply(obs.det.mat, 1, identical, x))})
      exp.hists <- apply(obs.det.hists, 1, FUN=calc.hist.prob, data=visit.frame, psi.table=psi.tab.new)
      
      visit.stats[j] <- sum(((obs.hists.counts-exp.hists)^2)/exp.hists) + (sum(obs.hists.counts)-sum(exp.hists))
    }
    grp.stats[i] <- sum(visit.stats)
  }
  return(test.stat=sum(grp.stats))
}

##Final function performs the goodness-of-fit permuation test
gof.test <-function(occ.form, det.form, model.fit, nsims, det.matrix, site.covs, group.covs, visit.covs){
  sim.dets <- simulate.det.hist(model.fit, nsims)
  gof.sims <- apply(sim.dets, 1, function(x){
    calc.gof.stat(fit.occ.model(occ.form, det.form, det.matrix, site.covs, 
                                group.covs, visit.covs,
                                mark.corr=model.fit$mark.corr, use.det.vec=TRUE, det.vec=x))
  })
  obs.stats <- calc.gof.stat(model.fit)
  pval1 <- mean(obs.stats < unlist(gof.sims))
  
  return(list(obs.stat = obs.stats,
              p.val = pval1))
}