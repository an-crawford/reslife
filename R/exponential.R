###########################
##Author: Andrew Crawford##
###########################

exp.rl <- function(fsoutput,x, p=.5, type = 'all') {
  if (fsoutput$ncovs == 0) {
    lambda = 1/exp(fsoutput$coefficients)
  }
  else {
    lambda = 1/exp(as.matrix(fsoutput$data$mml$rate) %*% as.numeric(fsoutput$coefficients));
  }
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

