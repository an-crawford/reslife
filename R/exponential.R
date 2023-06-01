data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ group, data = bc, dist = 'exp')
fsr$coefficients
fsr$res[1][1]
fsr


exp.rl <- function(fsoutput,x, p=.5, type = 'all') {
  lambda = 1/exp(as.matrix(fsoutput$data$mml$rate) %*% as.numeric(fsoutput$coefficients));
  mx = lambda
  sx = exp(-(1/lambda)*x)
  px = function(p){
    pc = (1-p)*sx
    px = qexp(pc, 1/lambda, lower.tail = FALSE)
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
exp.rl(fsr, 12, .5, 'mean')
