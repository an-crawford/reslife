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

weibull_mlr <- function(fsoutput,life, p=.5, type){

  if (fsoutput$dlist$name != "weibull.quiet") {
    print("wrong distribution")
  }

  else  {

    para_vect = fsoutput$coefficients
    shape_para = exp(para_vect[1])
    if (fsoutput$ncovs == 0){
      scale_para = exp(para_vect[2])
    }
    else{
      scale_para = exp(as.matrix(fsoutput$data$mml$scale) %*% as.numeric(para_vect[c(2:length(para_vect))]))
    }
    life_vec = rep(life, times = length(scale_para))
    zx = life^shape_para
    #S_zt = exp(-((zx/scale_para)^shape_para))
    S_t = exp(-((life/scale_para)^shape_para))
    #S_t = -scale*(life^shape)
    #mx <- as.numeric((((scale_para^shape_para)/shape_para)*gamma(1/shape_para)*S_zt)/(S_t))
    if (fsoutput$ncovs == 0){
      inte = function(x) {exp(-((x/scale_para)^shape_para))}
      v = integrate(inte, lower = life, upper= Inf)
      mx = as.numeric(v$value/S_t)
    }
    else{
      mx = vector('numeric', length = length(scale_para))
      for (i in 1:length(scale_para)){
        integ = function(x) {exp(-((x/scale_para[i])^shape_para))}
        v = integrate(integ, lower = life, upper= Inf)
        mx[i] = v$value/S_t[i]
      }
    }

    px = function(p){
      pc = (1-p)*S_t
      px = as.numeric(qweibull(pc, shape = shape_para, scale = scale_para, lower.tail = FALSE) - life)
    }

    if (type=='mean'){
      return(list('mean'= mx))
    }
    if (type== 'median'){
      return(list('median'= px(.5)))
    }
    if (type == 'percentile'){
      return(list('percentile'= px(p)))
    }
    if (type == 'all'){
      return(list("mean" = mx, 'median' = px(.5), 'percentile' = px(p)))
    }
    else{
      return('invalid type')
    }

  }

}

