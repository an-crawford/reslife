
data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ group, data = bc, dist = 'gamma')
length(fsr$coefficients)
s= fsr$res[2][1]
s =100

gamma.rl = function(fsr, x, p=.5, type = 'all'){
  a = fsr$res[1][1]
  lambda = 1/(fsr$res[2][1])
  sx = pgamma(x, shape = a, scale = lambda, lower.tail = FALSE)
  mx = ((x^a)*exp(-x/lambda))/((lambda^(a-1))*gamma(a)*sx)+(lambda*a)
  px = function(p){
    pc = (1-p)*sx
    px = qgamma(pc, shape = a, scale = lambda, lower.tail = FALSE)
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
gamma.rl(fsr, 10, .99, 'all')
