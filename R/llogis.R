###########################
##Author: Andrew Crawford##
###########################

llogis.rl = function(fsroutput, x, p=.5, type = 'all'){
  a = exp(fsroutput$coefficients[1])
  if (fsroutput$ncovs == 0) {
    lambda = as.numeric(exp(fsroutput$coefficients[2]))
  }
  else{
    s = fsr$coefficients
    #s[2] = exp(s[2])
    lambda = exp(as.matrix(fsroutput$data$mml$scale) %*% as.numeric(s[-1]))
  }
  d = (x/lambda)^a
  sx = 1/(1+d)
  zx = (d)/(1+(d))
  shapez = 1- (1/a)
  scalez = 1/a
  sz = pbeta(zx, shape1 = shapez, shape2 = scalez, lower.tail = FALSE)
  #mz = as.numeric((lambda/a)*gamma(1-(1/a))*gamma(1/a)*sz*(1+d))
  mx = as.numeric((lambda*(pi/a)/sin(pi/a))*(1+d)*(zipfR::Rbeta((1-sx), scalez, shapez, lower = FALSE)))
  px = function(p){
    pc = (1-p)*sx
    px = as.numeric(qllogis(pc, shape = a, scale = (lambda), lower.tail = FALSE) - x)
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
library(flexsurv)
fsr = flexsurvreg(Surv(recyrs, censrec)~ group, data = bc, dist = 'llogis')
llogis.rl(fsr, 7, .99, 'mean')
fsr$dlist$name
