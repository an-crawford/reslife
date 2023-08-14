###########################
##Author: Andrew Crawford##
###########################


lnorm_rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){
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
  sigma = exp(fsroutput$coefficients[2])
  if (length(newdata) == 0){
    if (fsroutput$ncovs == 0) {
      mu = fsroutput$coefficients[1]
    }
    else{
      mu = (as.matrix(fsroutput$data$mml$meanlog) %*% as.numeric(fsroutput$coefficients[-2]))
    }
  }
  else{
    if (fsroutput$ncovs == 0) {
      mu = fsroutput$coefficients[1]
    }
    else{
      X<-model.matrix( ~ ., data = newdata)
      s = fsroutput$coefficients
      sa = s[1]
      sb = s[-c(1,2)]
      sb = sb[colnames(X)]
      sb = sb[!is.na(sb)]
      sc = append(sa,sb)
      if (length(sc) != ncol(X)){
        stop('Incorrect Level Entered')
      }
      mu = (as.matrix(X) %*% as.numeric(sc))
    }
  }
  #mu <- exp(meanlog + (1/2)*sdlog^2)
  #sigma <- exp(meanlog + (1/2)*sdlog^2)*sqrt(exp(sdlog^2) - 1)
  sx = pnorm(log(x), mean = mu, sd = sigma, lower.tail=FALSE)
  mx = as.numeric(exp(mu + (sigma^2)/2)*(pnorm(log(x), mean = (mu + sigma^2), sd = sigma, lower.tail = FALSE))/sx - x)
  px = function(p){
    pc = (1-p)*sx
    px = as.numeric(qlnorm(pc, meanlog = mu, sdlog = sigma, lower.tail = FALSE)-x)
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
