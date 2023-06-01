#install.packages('flexsurv')
#install.packages('survival')
library(flexsurv)

#input and output
#input
#flexresurv(fsoutput, <- object
#           life), <- float > 0

#output numerics

#inside of main function
#fs <- flexsurvreg(Surv(response_var, censor_indcator) ~ variable_vector, 
#                 data = data_name,
#                 dist = dist_name)

#if (type == 'mean') { 
# result = MeanResLifeFunc (fs, dist_name, life)
# } else if (type == 'median') {
#   result = MedianResLifeFunc (fs, dist_name, life)
# } else if  (type =='percentile') {
#   result = PercResLifeFunc (fs, dist_name, life, percentile)
# } else {
#   'report a type error'
# } 
#return result


### Written by Zekai Wang ###

weibull_mlr <- function(fsoutput,life){
  
  if (fsoutput$dlist$name != "weibull.quiet") {
    print("wrong distribution")
  }
  
  else  {
    
    para_vect = fsoutput$coefficients
    shape_para = para_vect[1]
    scale_para = para_vect[2]
    zx = life^shape_para
    S_zt = exp(-((zx/scale_para)^shape_para))
    S_t = exp(-((life/scale_para)^shape_para))
    #S_t = -scale*(life^shape)
    
    if (fsoutput$ncovs == 0){
      
      mx <- as.numeric((((scale_para^shape_para)/shape_para)*gamma(1/shape_para)*S_zt)/(S_t))
      
    }
    else{
      
      scale_para = exp(as.matrix(fsoutput$data$mml$scale) %*% as.numeric(para_vect[c(2:length(para_vect))]))
      mx <- as.numeric((((scale_para^shape_para)/shape_para)*gamma(1/shape_para)*S_zt)/(S_t))
      
    }
    
    
    return (mx)
    
  }
  
}

#test case


fse <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc, dist = "exp")
fsw <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc, dist = "weibull")
fsw2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc, dist = "weibull")

weibull_mlr(fs_s, 1)
