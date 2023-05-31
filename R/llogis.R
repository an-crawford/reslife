data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ 1, data = bc, dist = 'llogis')
length(fsr$coefficients)
fsr$res[1][1]
fsr



llogis.rl = function(fsr, x, p=.5, type = 'all'){
  a = fsr$res[1][1]
  lambda = fsr$res[2][1]
  sx = (1+(x/lambda)^a)^-1
  zx = ((x/lambda)^a)/(1+(x/lambda)^a)
  shapez = 1- (1/a)
  scalez = 1/a
  sz = (1+(zx/scalez)^shapez)^-1
  mx = x+(lambda/a)*gamma(1-(1-a))*gamma(1/a)*sz*(1/sx)
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
llogis.rl(fsr, 14, .875, 'all')
