###########################
##Author: Andrew Crawford##
###########################

exp_rl <- function(fsoutput,x, p=.5, type = 'all', newdata = data.frame()) {
  if (length(newdata)!=0){
    if (length(newdata) == 1){
      if(fsoutput$covdata$covnames != colnames(newdata)){
        stop('Wrong columns in inputted data')
      }
    }
    else{
      names = fsoutput$covdata$covnames
      newdata= newdata[,c(names)]
    }
  }


  if (length(newdata) == 0){

    if (fsoutput$ncovs == 0) {
      lambda = 1/exp(fsoutput$coefficients)
    }
    else {
      lambda = 1/exp(as.matrix(fsoutput$data$mml$rate) %*% as.numeric(fsoutput$coefficients))
    }
  }
  else{
    if (fsoutput$ncovs == 0) {
      lambda = 1/exp(fsoutput$coefficients)
    }
    else{
      X<-model.matrix( ~ ., data = newdata)
      s = fsoutput$coefficients
      sa = s[1]
      sb = s[-1]
      sb = sb[colnames(X)]
      sb = sb[!is.na(sb)]
      sc = append(sa,sb)
      if (length(sc) != ncol(X)){
        stop('Incorrect Level Entered')
      }
      lambda = 1/exp(as.matrix(X) %*% as.numeric(sc))
    }
  }
  mx = lambda
  sx = exp(-(1/lambda)*x)
  px = function(p){
    pc = (1-p)*sx
    px = qexp(pc, 1/lambda, lower.tail = FALSE)
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

