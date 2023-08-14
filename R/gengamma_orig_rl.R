#######################
#### GenGamma Orig ####
#######################

####Written by Zekai Wang ####


upper_incomplete_gamma <- function(x,a) {
  #return (incgam(x,a))
  return (gamma(a) * pgamma(x, a, 1, lower.tail = FALSE))
}


gengamma_orig_rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){
  if (length(newdata)!=0){
    if (length(newdata) == 1){
      if(fsroutput$covdata$covnames != colnames(newdata)){
        stop('Wrong columns in inputted data')
      }
    }
    else{
      names = fsroutput$covdata$covnames
      newdata= newdata[,c(names)]
      #print(newdata)
    }
  }
  b = as.numeric(exp(fsroutput$coefficients[1]))
  k = as.numeric(exp(fsroutput$coefficients[3]))
  if (length(newdata) == 0){
    if (fsroutput$ncovs == 0) {
      a = as.numeric(exp(fsroutput$coefficients[2]))
    }
    else{
      s = fsroutput$coefficients
      a = exp(as.matrix(fsroutput$data$mml$scale) %*% as.numeric(s[-c(1,3)]))
    }
  }
  else{
    if (fsroutput$ncovs == 0) {
      a = as.numeric(exp(fsroutput$coefficients[2]))
    }
    else{
      X<-model.matrix( ~ ., data = newdata)
      s = fsroutput$coefficients
      sa = s[2]
      sb = s[-c(1,2,3)]
      sb = sb[colnames(X)]
      sb = sb[!is.na(sb)]
      sc = append(sa,sb)
      if (length(sc) != ncol(X)){
        stop('Incorrect Level Entered')
      }
      a = exp(as.matrix(X) %*% as.numeric(sc))
    }
  }

  sx = pgengamma.orig(x, shape = b, scale = a, k = k, lower.tail = FALSE)
  mx = as.numeric(a * upper_incomplete_gamma((x/a)^b, k + 1/b)/upper_incomplete_gamma((x/a)^b, k) - x)
  px = function(p){
    pc = (1-p)*sx
    px = as.numeric(qgengamma.orig(pc, shape = b, scale = a, k = k, lower.tail = FALSE) - x)
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
