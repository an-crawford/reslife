###Written by Ka Lok Lee

###Load flexsurv
library(flexsurv)

###Check dataset bc
head(bc)

###Run Weibull
fsw<-flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc,dist = "weibull")
fsw$coefficients

###Run Exponential
fse<-flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc,dist = "exp")
fse$coefficients

###Exponential MRL function
exp_mrl <- function(fsoutput,life) {
  
  if (fsoutput$dlist$name != 'exp') {
    print("error")
  }
  
  else  {

  mean_res_life = 1/exp(as.matrix(fsoutput$data$mml$rate) %*% as.numeric(fsoutput$coefficients));
  
  return (mean_res_life)
  
  }
}

###Calculate MRL
exp_mrl(fse,1)

###Calculate using the wrong object
exp_mrl(fsw,1)

