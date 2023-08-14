
###########################
##Author: Andrew Crawford##
###########################

gamma_rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){
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
        lambda = 1/as.numeric(exp(fsroutput$coefficients[2]))
      }
      else{
        s = fsroutput$coefficients
        lambda = 1/exp(as.matrix(fsroutput$data$mml$rate) %*% as.numeric(s[-1]))
      }
    }
    else{
      if (fsroutput$ncovs == 0) {
        lambda = 1/as.numeric(exp(fsroutput$coefficients[2]))
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
        lambda = 1/exp(as.matrix(X) %*% as.numeric(sc))
      }
    }

  #lambda = 1/(exp(as.matrix(fsroutput$data$mml$rate) %*% as.numeric(fsroutput$coefficients[-1])))
  sx = pgamma(x, shape = a, scale = lambda, lower.tail = FALSE)
  mx = as.numeric(((x^a)*exp(-x/lambda))/((lambda^(a-1))*gamma(a)*sx)+(lambda*a) - x)
  px = function(p){
    pc = (1-p)*sx
    px = as.numeric(qgamma(pc, shape = a, scale = lambda, lower.tail = FALSE) - x)
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

