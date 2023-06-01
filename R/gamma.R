
data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ 1, data = bc, dist = 'gamma')
fsr$coefficients[1]
s= fsr$res[2][1]
s =100
fsr$dlist$name
lambda = 1/(exp(as.matrix(fsr$data$mml$rate) %*% as.numeric(fsr$coefficients[-1])))


gamma.rl = function(fsroutput, x, p=.5, type = 'all'){
  a = exp(fsroutput$coefficients[1])
  lambda = 1/(exp(as.matrix(fsroutput$data$mml$rate) %*% as.numeric(fsroutput$coefficients[-1])))
  sx = pgamma(x, shape = a, scale = lambda, lower.tail = FALSE)
  mx = ((x^a)*exp(-x/lambda))/((lambda^(a-1))*gamma(a)*sx)+(lambda*a)
  px = function(p){
    pc = (1-p)*sx
    px = qgamma(pc, shape = a, scale = lambda, lower.tail = FALSE)
  }
  if (type=='mean'){
    colnames(mx) = 'mean'
    return(list('mean'= unique(cbind(fsroutput$data$mml$rate, mx))))
  }
  if (type== 'median'){
    med = px(.5)
    colnames(med) = 'median'
    return(list('median'= unique(cbind(fsroutput$data$mml$rate, med))))
  }
  if (type == 'percentile'){
    per = px(p)
    #colnames(per) = cat(p, 'percentile')
    colnames(per) = sprintf('%s percentile', p)
    return(list('percentile'= unique(cbind(fsroutput$data$mml$rate, per))))
  }
  if (type == 'all'){
    colnames(mx) = 'mean'
    med = px(.5)
    colnames(med) = 'median'
    per = px(p)
    colnames(per) = sprintf('%s percentile',p)
    return(list('all'= unique(cbind(fsroutput$data$mml$rate, mx, med, per))))
  }
  else{
    return('invalid type')
  }


}
gamma.rl(fsr, 10, .99, 'mean')

