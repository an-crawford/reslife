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
    s[2] = exp(s[2])
    lambda = (as.matrix(fsroutput$data$mml$scale) %*% as.numeric(s[-1]))
  }
  sx = (1+(x/lambda)^a)^-1
  zx = ((x/lambda)^a)/(1+(x/lambda)^a)
  shapez = 1- (1/a)
  scalez = 1/a
  sz = (1+(zx/scalez)^shapez)^-1
  mx = x+(lambda/a)*gamma(1-(1/a))*gamma(1/a)*sz*(1/sx)
  px = function(p){
    pc = (1-p)*sx
    px = qllogis(pc, shape = a, scale = (lambda), lower.tail = FALSE)
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

