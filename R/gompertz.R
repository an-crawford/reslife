data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ group, data = bc, dist = 'gompertz')
fsr$coefficients[1]


gompertz.rl = function(fsroutput, x, p=.5, type = 'all'){
  a = fsroutput$coefficients[1]
  lambda = (exp(as.matrix(fsroutput$data$mml$rate) %*% as.numeric(fsroutput$coefficients[-1])))
  sx = exp(((lambda)/a)*(1-exp(a*x)))
  zx = ((lambda)/a)*exp(a*x)
  mx = x+exp(zx)*(1/a)*Vectorize(incgam)(zx, 0)
  px = function(p){
    pc = (1-p)*sx
    px = qgompertz(pc, shape = a, rate = (lambda), lower.tail = FALSE)
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
gompertz.rl(fsr, 14, .75, 'percentile')

