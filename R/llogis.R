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
  sx = (1+(x/lambda)^a)^-1
  #zx = ((x/lambda)^a)/(1+(x/lambda)^a)
  #shapez = 1- (1/a)
  #scalez = 1/a
  #sz = (1+(zx/scalez)^shapez)^-1
  #mx = as.numeric((lambda/a)*gamma(1-(1/a))*gamma(1/a)*sz*(1/sx))
  if (fsroutput$ncovs == 0){
    inte = function(t) {(1+(t/lambda)^a)^-1}
    v = integrate(inte, lower = x, upper= Inf)
    mx = as.numeric(v$value/sx)
  }
  else{
    mx = vector('numeric', length = length(lambda))
    for (i in 1:length(lambda)){
      integ = function(t) {(1+(t/lambda[i])^a)^-1}
      v = integrate(integ, lower = x, upper= Inf)
      mx[i] = v$value/sx[i]
    }
  }
  px = function(p){
    pc = (1-p)*sx
    px = as.numeric(qllogis(pc, shape = a, scale = (lambda), lower.tail = FALSE) - x)
  }
  if (type=='mean'){
    return(list('mean'= mx, 'inte' = mz))
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
fsr = flexsurvreg(Surv(recyrs, censrec)~ group, data = bc, dist = 'llogis')
residLife(fsr, 2, p = .95, type = 'all')
