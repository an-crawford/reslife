#install.packages('flexsurv')
#install.packages('survival')
library(survival)
library(flexsurv)

#input and output
#input
#flexresurv(response_var, <- string
#           censor_indcator, <- string
#           variable_vector, <- vector of strings, 1 for non-variable
#           data_name, <- string, dataframe name
#           dist_name, <- string
#           type <- string, what type of residual lifetime function
#           life), <- float > 0

#output mlr(life)

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

fse <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc, dist = "exp")
fsw <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc, dist = "weibull")

weibull_mlr <- function(fsoutput,life){
  
  if (fsoutput$dlist$name != "weibull.quiet") {
    print("error")
  }
  
  else  {
    para_vect = fsoutput$coefficients
    shape = para_vect[1]
    scale = para_vect[2]
    zx = life^shape
    S_zt = exp(-(zx/scale)^shape)
    S_t = exp(-(life/scale)^shape)
    scale_para = exp(as.matrix(fsoutput$data$mml$scale) %*% as.numeric(para_vect[c(2:length(para_vect)-1)]))
    mx = as.numeric(((scale_para/shape)*gamma(1/shape)*S_zt)/(S_t))
    
    return (mx)
    
  }
  
}

weibull_mlr(fsw, 10)
