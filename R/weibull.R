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

weibull_mlr <- function(fsoutput,life, p=.5, type = 'all', newdata = data.frame()){

  if (fsoutput$dlist$name != "weibull.quiet") {
    print("wrong distribution")
  }

  else  {
    if (length(newdata)!=0){
      if (length(newdata) == 1){
        stopifnot(fsoutput$covdata$covnames == colnames(newdata))
      }
      else{
        names = fsoutput$covdata$covnames
        newdata= newdata[,c(names)]
        #print(newdata)
      }
    }
    para_vect = fsoutput$coefficients
    shape_para = exp(para_vect[1])
    if (length(newdata) == 0){
      if (fsoutput$ncovs == 0){
        scale_para = exp(para_vect[2])
      }
      else{
        scale_para = exp(as.matrix(fsoutput$data$mml$scale) %*% as.numeric(para_vect[c(2:length(para_vect))]))
      }
    }
    else{
      if (fsoutput$ncovs == 0){
        scale_para = exp(para_vect[2])
      }
      else{
        X<-model.matrix( ~ ., data = newdata)
        s = fsoutput$coefficients
        sa = s[2]
        sb = s[-c(1,2)]
        sb = sb[colnames(X)]
        sb = sb[!is.na(sb)]
        sc = append(sa,sb)
        if (length(sc) != ncol(X)){
          print('Incorrect Level Entered')
          error = 1
          stopifnot(error = 0)
        }
        scale_para = exp(as.matrix(X) %*% as.numeric(sc))
      }
    }
    S_t = exp(-((life/scale_para)^shape_para))
    mx = as.numeric((scale_para/shape_para)*exp((life/scale_para)^shape_para)*Vectorize(pracma::incgam)((life/scale_para)^shape_para, 1/shape_para))

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

