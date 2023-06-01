data("bc")
fsr = flexsurvreg(Surv(recyrs, censrec)~ 1, data = bc, dist = 'lnorm')
fsr$coefficients
fsr

lnorm.rl = function(fsroutput, x, p=.5, type = 'all'){
  sigma = fsroutput$coefficients[2]
  mu = (exp(as.matrix(fsroutput$data$mml$meanlog) %*% as.numeric(fsroutput$coefficients[-2])))
  sx = pnorm(x, mean = exp(mu), sd = exp(sigma), lower.tail=FALSE)
  mx = log(exp(mu + (sigma^2)/2)*(pnorm(x, mean = exp(mu + sigma^2), sd = exp(sigma), lower.tail = FALSE))/sx)
  px = function(p){
    pc = (1-p)*sx
    px = log(qnorm(pc, mean = exp(mu), sd = exp(sigma), lower.tail = FALSE))
  }
  if (type=='mean'){
    colnames(mx) = 'mean'
    return(list('mean'= unique(cbind(fsroutput$data$mml$meanlog, mx))))
  }
  if (type== 'median'){
    med = px(.5)
    colnames(med) = 'median'
    return(list('median'= unique(cbind(fsroutput$data$mml$meanlog, med))))
  }
  if (type == 'percentile'){
    per = px(p)
    #colnames(per) = cat(p, 'percentile')
    colnames(per) = sprintf('%s percentile', p)
    return(list('percentile'= unique(cbind(fsroutput$data$mml$meanlog, per))))
  }
  if (type == 'all'){
    colnames(mx) = 'mean'
    med = px(.5)
    colnames(med) = 'median'
    per = px(p)
    colnames(per) = sprintf('%s percentile',p)
    return(list('all'= unique(cbind(fsroutput$data$mml$meanlog, mx, med, per))))
  }
  else{
    return('invalid type')
  }
}
lnorm.rl(fsr, 10, .99, 'mean')


sigma = fsr$coefficients[2]
mu = (exp(as.matrix(fsr$data$mml$meanlog) %*% as.numeric(fsr$coefficients[-2])))
sx = pnorm(12, mean = exp(mu), sd = exp(sigma), lower.tail=FALSE)
mu
sigma
exsx
log(100)
