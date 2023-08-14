###########################
##Author: Andrew Crawford##
###########################

llogis_rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){
  regbeta = function (x, a, b, lower = TRUE, log = !missing(base), base = exp(1))
  {
    if (log) {
      if (missing(base)) {
        pbeta(x, shape1 = a, shape2 = b, lower.tail = lower,
              log.p = TRUE)
      }
      else {
        pbeta(x, shape1 = a, shape2 = b, lower.tail = lower,
              log.p = TRUE)/log(base)
      }
    }
    else {
      pbeta(x, shape1 = a, shape2 = b, lower.tail = lower,
            log.p = FALSE)
    }
  }
  if (length(newdata)!=0){
    if (length(newdata) == 1){
      if(fsroutput$covdata$covnames != colnames(newdata)){
        stop('Wrong columns in inputted data')
      }
    }
    else{
      names = fsroutput$covdata$covnames
      newdata= newdata[,c(names)]
    }
  }
  a = exp(fsroutput$coefficients[1])
  if (length(newdata) == 0){
    if (fsroutput$ncovs == 0) {
      lambda = as.numeric(exp(fsroutput$coefficients[2]))
    }
    else{
      s = fsroutput$coefficients
      #s[2] = exp(s[2])
      lambda = exp(as.matrix(fsroutput$data$mml$scale) %*% as.numeric(s[-1]))
    }
  }
  else{
    if (fsroutput$ncovs == 0) {
      lambda = as.numeric(exp(fsroutput$coefficients[2]))
    }
    else{
      X<-model.matrix( ~ ., data = newdata)
      s = fsroutput$coefficients
      sa = s[2]
      sb = s[-c(1,2)]
      sb = sb[colnames(X)]
      sb = sb[!is.na(sb)]
      sc = append(sa,sb)
      if (length(sc) != ncol(X)){
        stop('Incorrect Level Entered')
      }
      lambda = exp(as.matrix(X) %*% as.numeric(sc))
    }
  }
  d = (x/lambda)^a
  sx = 1/(1+d)
  zx = (d)/(1+(d))
  shapez = 1- (1/a)
  scalez = 1/a
  sz = pbeta(zx, shape1 = shapez, shape2 = scalez, lower.tail = FALSE)
  #mz = as.numeric((lambda/a)*gamma(1-(1/a))*gamma(1/a)*sz*(1+d))
  mx = as.numeric((lambda*(pi/a)/sin(pi/a))*(1+d)*(regbeta((1-sx), scalez, shapez, lower = FALSE)))
  px = function(p){
    pc = (1-p)*sx
    px = as.numeric(qllogis(pc, shape = a, scale = (lambda), lower.tail = FALSE) - x)
  }
  if (type=='mean'){
    return(c(mx))
  }
  if (type== 'median'){
    return(c(px(.5)))
  }
  if (type == 'percentile'){
    return(c(px(p)))
  }
  if (type == 'all'){
    return(data.frame(mean = mx, median = px(.5), percentile = px(p)))
  }
  else{
    return('invalid type')
  }
}
