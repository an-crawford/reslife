
###########################
##Author: Andrew Crawford##
###########################

gamma.rl = function(fsroutput, x, p=.5, type = 'all', newdata = c()){
  if (fsroutput$covdata$isfac == TRUE){
    lvl = fsroutput$covdata$xlev$group
    #return(newdata)
    stopifnot(newdata %in% lvl)
  }
  a = exp(fsroutput$coefficients[1])
  if (length(newdata) == 0){
    if (fsroutput$ncovs == 0) {
      lambda = 1/as.numeric(exp(fsroutput$coefficients[2]))
    }
    else{
      s = fsroutput$coefficients
      lambda = 1/exp(as.matrix(fsroutput$data$mml$rate) %*% as.numeric(s[-1]))
    }
  }
  else{
    if (fsroutput$ncovs == 0) {
      lambda = 1/as.numeric(exp(fsroutput$coefficients[2]))
    }
    else{
      s = fsroutput$coefficients
      b = model.matrix(~0+ newdata)
      lambda = 1/exp(as.matrix(b) %*% as.numeric(s[-1]))
    }
  }
  #lambda = 1/(exp(as.matrix(fsroutput$data$mml$rate) %*% as.numeric(fsroutput$coefficients[-1])))
  sx = pgamma(x, shape = a, scale = lambda, lower.tail = FALSE)
  mx = as.numeric(((x^a)*exp(-x/lambda))/((lambda^(a-1))*gamma(a)*sx)+(lambda*a) - x)
  px = function(p){
    pc = (1-p)*sx
    px = as.numeric(qgamma(pc, shape = a, scale = lambda, lower.tail = FALSE) - x)
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


fsr = flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc,dist = "gamma")
newdat = c(group = 'Good', 'Poor')
#nd = as.data.frame(newdat)
#nd$V3 = as.factor(nd$V3)
gamma.rl(fsr, 4, type = 'mean', newdata = newdat)
#predict(fsr, newdata = as.data.frame(newdat))
