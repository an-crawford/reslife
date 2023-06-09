
###########################
##Author: Andrew Crawford##
###########################

gamma.rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){
  if (length(newdata)!=0){
    if (length(newdata) == 1){
      stopifnot(fsroutput$covdata$covnames == colnames(newdata))
    }
    else{
      names = fsroutput$covdata$covnames
      newdata= newdata[,c(names)]
      print(newdata)
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
        print('Incorrect Level Entered')
        error = 1
        stopifnot(error = 0)
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

group = c("Medium", 'Good', "Poor")
age = c(43, 35, 39)
newdata = data.frame(age, group)
newd = data.frame(group)

fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gamma")
fsr = flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc,dist = "gamma")
newbc <- bc
newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE),sd = 5)
fsr2 =  flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc,
                    dist="gamma")


gamma.rl(fitg, 4, type = 'mean', newdata = newd)
gamma.rl(fsr2, 4, type = 'mean', newdata = newdata)
#predict(fsr, newdata = as.data.frame(newdat))

fsr2
name = fsr2$covdata$covnames
newdata2= data.frame(newdata[,c(name)])
#colnames(newdata2) = c(name)
newdata2
name
