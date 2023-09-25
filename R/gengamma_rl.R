####################
##### GenGamma #####
####################

############################
#### Author: Zekai Wang 
####Updated by Ka Lok Lee Sep 2023
############################

#When Q > 0, the parameters relationship between gengamma and gengamma.orig is:
#dgengamma.orig(x, shape=abs(Q)/(sigma), scale=exp(mu - (log(Q^(-2))*sigma)/(abs(Q))), k=1/(Q^2)) == dgengamma(x, mu=mu, sigma=sigma, Q=abs(Q))
# ------------------------------------------
#dgengamma.orig(x, shape=b, scale=a, k) == dgengamma(x, mu=log(a) + log(k)/b, sigma=1/(b*sqrt(k)), Q=1/sqrt(k))
######### example #######
#library('flexsurv')
#fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc ,dist = "genGamma")
#fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc ,dist = "genGamma.orig")

#mu <- as.numeric(fs1$coefficients[1])
#sigma <- exp(as.numeric(fs1$coefficients[2]))
#Q <- as.numeric(fs1$coefficients[3])
#x = 4
#dgengamma.orig(x, shape=abs(Q)/(sigma), scale=exp(mu - (log(Q^(-2))*sigma)/(abs(Q))), k=1/(Q^2))
#dgengamma(x, mu=mu, sigma=sigma, Q=abs(Q))

#b <- exp(as.numeric(fs2$coefficients[1]))
#a <- exp(as.numeric(fs2$coefficients[2]))
#k <- exp(as.numeric(fs2$coefficients[3]))

#dgengamma.orig(x, shape=b, scale=a, k)
#dgengamma(x, mu=log(a) + log(k)/b, sigma=1/(b*sqrt(k)), Q=1/sqrt(k))
# ---------------------------------------------
# ---------------------------------------------


upper_incomplete_gamma <- function(a,x) {
  
  return (gamma(a) * pgamma(x, a, 1, lower.tail = FALSE))
  
}

lower_incomplete_gamma <- function(a,x) {
  
  return (gamma(a) - upper_incomplete_gamma(a,x))
  
}

gengamma_rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){
  
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
  
  sigma = exp(as.numeric(fsroutput$coefficients[2]))
  Q = as.numeric(fsroutput$coefficients[3])
  
if ( Q > 0 ) {  ##### Transform the GenGamma parameters into GenGamma.orig and use that MRL
  
  k = 1 / (Q^2)
  b = abs(Q)/(sigma)
  
  if (length(newdata) == 0){
    
    if (fsroutput$ncovs == 0) {
      
      mu = as.numeric(fsroutput$coefficients[1])
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))
      
    }
    
    else{
      
      s = fsroutput$coefficients
      
      mu = as.matrix(fsroutput$data$mml$mu) %*% as.numeric(s[-c(2,3)])
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))
    }
    
  }
  
  else{
    
    if (fsroutput$ncovs == 0) {
      
      mu = as.numeric(fsroutput$coefficients[1])
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))
      
      
    }
    
    else{
      
      X<-model.matrix( ~ ., data = newdata)
      
      s = fsroutput$coefficients
      
      sa = s[1]
      
      sb = s[-c(1,2,3)]
      
      sb = sb[colnames(X)]
      
      sb = sb[!is.na(sb)]
      
      sc = append(sa,sb)
      
      if (length(sc) != ncol(X)){
        
        stop('Incorrect Level Entered')
        
        
      }
      
      mu = as.matrix(X) %*% as.numeric(sc)
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))
      
    }
    
  }
  
  
  mx = as.numeric(a * upper_incomplete_gamma(k + 1/b, (x/a)^b) / upper_incomplete_gamma(k, (x/a)^b) - x)

}
  
  
if (Q < 0) {  
  

  if (length(newdata) == 0){
    
    if (fsroutput$ncovs == 0) {
      
      mu = as.numeric(fsroutput$coefficients[1])
      
    }
    
    else{
      
      mu = as.matrix(fsroutput$data$mml$mu) %*% as.numeric(fsroutput$coefficients[-c(2,3)])

    }
    
  }
  
  else{
    
    if (fsroutput$ncovs == 0) {
      
      mu = as.numeric(fsroutput$coefficients[1])
      
    }
    
    else{
      
      X<-model.matrix( ~ ., data = newdata)
      
      s = fsroutput$coefficients
      
      sa = s[1]
      
      sb = s[-c(1,2,3)]
      
      sb = sb[colnames(X)]
      
      sb = sb[!is.na(sb)]
      
      sc = append(sa,sb)
      
      if (length(sc) != ncol(X)){
        
        stop('Incorrect Level Entered')
        
        
      }
      
      mu = as.matrix(X) %*% as.numeric(sc)
        
    }
    
  }
  
  Cx = Q^(-2)*(exp(-mu)*x)^(Q/sigma)
  
  mx = as.numeric( ( exp(mu)*(Q^2)^(sigma/Q)*lower_incomplete_gamma(Q^(-2)+sigma/Q,Cx)-
                       x*lower_incomplete_gamma(Q^(-2),Cx) ) / lower_incomplete_gamma(Q^(-2),Cx) )
  
}
  
  sx = pgengamma(x,
                 
                 mu = mu,
                 
                 sigma = sigma,
                 
                 Q = Q,
                 
                 lower.tail = FALSE)
  
  
  
  px = function(p){
    
    pc = (1-p)*sx
    
    px = as.numeric(qgengamma(pc,
                              
                              mu = mu,
                              
                              sigma = sigma,
                              
                              Q = Q,
                              
                              lower.tail = FALSE) - x)
    
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




# Accuracy Test #
### Q>0
#fs1 <- flexsurvreg(Surv(years, status) ~ 1, data = bosms3 ,dist = "genGamma")
#mrl <- gengamma_rl(fs1, 4, p=.5, type ='mean')
#mu <- as.numeric(fs1$coefficients[1])
#sigma <- exp(as.numeric(fs1$coefficients[2]))
#Q <- as.numeric(fs1$coefficients[3]) ### Q>0
##f1 <- function(x) {
#return(pgengamma(x,mu = mu,sigma = sigma,Q = Q, lower.tail = FALSE))
#}
#ni <- integrate(f1, 4, Inf)$value/pgengamma(4, mu = mu,sigma = sigma,Q = Q, lower.tail = FALSE)
#ni
#mrl

### Q<0
#newbc <- bc
#newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE),sd = 5)
#fsr2 <-  flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc, dist = 'gengamma')
#mrl <- gengamma_rl(fsr2, 4, p=.5, type ='mean')
#mu = as.matrix(fsr2$data$mml$mu) %*% as.numeric(fsr2$coefficients[-c(2,3)])
#sigma <- exp(as.numeric(fsr2$coefficients[2]))
#Q <- as.numeric(fsr2$coefficients[3]) ### Q<0
#f1 <- function(x) {
# return(pgengamma(x,mu = mu[1],sigma = sigma,Q = Q, lower.tail = FALSE))
#}
#ni <- integrate(f1, 4, Inf)$value/pgengamma(4, mu = mu[1],sigma = sigma,Q = Q, lower.tail = FALSE)
#ni
#mrl[1]
