####################

##### GenF.orig #####

####################


############################

#### Author: Zekai Wang ####

#### Updated by: Ka Lok Lee Sep 2023 ###

############################



# parameters relationship between Cox’s parameterization and genF.orig

# s1 = m1, s2 = m2, sigma = sigma, mu = beta


incomplete_beta <- function(x,a,b) {
  return (pbeta(x,a,b)*beta(a,b));
}

Gauss2F1b <- function(a,b,c,x){
  ifelse(x>=0 & x<1, hyperg_2F1(a,b,c,x),hyperg_2F1(a,c-b,c,1-1/(1-x))/(1-x)^a)
}


genF_orig_rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){
  
  
  
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
  
  m1 = exp(as.numeric(fsroutput$coefficients[3]))
  
  m2 = exp(as.numeric(fsroutput$coefficients[4]))
  
  
  
  
  
  if (length(newdata) == 0){
    
    
    
    if (fsroutput$ncovs == 0) {
      
      
      
      mu = as.numeric(fsroutput$coefficients[1])
      
      
      
    }
    
    
    
    else{
      
      
      
      s = fsroutput$coefficients
      
      mu = as.matrix(fsroutput$data$mml$mu) %*% as.numeric(s[-c(2,3,4)])
      
      
      
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
      
      
      
      sb = s[-c(1,2,3,4)]
      
      
      
      sb = sb[colnames(X)]
      
      
      
      sb = sb[!is.na(sb)]
      
      
      
      sc = append(sa,sb)
      
      
      
      if (length(sc) != ncol(X)){
        
        
        
        stop('Incorrect Level Entered')
        
        
        
        
        
        
        
      }
      
      
      
      mu = as.matrix(X) %*% as.numeric(sc)
      
      
      
    }
    
    
    
  }
  
  
  sx = pgenf.orig(x,
                  
                  
                  mu = mu,
                  
                  
                  sigma = sigma,
                  
                  
                  s1 = m1,
                  
                  
                  s2 = m2,
                  
                  
                  lower.tail = FALSE)
  
  
  C = exp(-mu/sigma)*(m1/m2)*(x^(1/sigma))
  
  part1 = (exp(mu)*((m2/m1)^sigma)*C^(-m2+sigma))/(beta(m1,m2))
  
  part2 = gamma(m2-sigma)/gamma(m2-sigma+1)
  
  part3 = Gauss2F1b(m1+m2, m2-sigma, m2-sigma+1, -1/C)
  
  num_integral = part1 * part2 * part3
  
###########Closed Form for raw mean#############################
  
  Genf_mean_raw <- function(mu, sigma, m1, m2) {

   if ( m2 < sigma ){
       message("Estimated parameters produce an undefined mean.")
   }
    return(exp(mu)*(m2/m1)^sigma*beta(m1+sigma,m2-sigma)/beta(m1,m2))
    
  }
  
  
  mx <- ifelse(C == 0, Genf_mean_raw(mu,sigma, m1, m2), num_integral/sx - x)
  
  px = function(p){
    
    
    
    pc = (1-p)*sx
    
    
    
    px = as.numeric(qgenf.orig(pc,
                               
                               
                               
                               mu = mu,
                               
                               
                               
                               sigma = sigma,
                               
                               
                               
                               s1 = m1,
                               
                               
                               
                               s2 = m2,
                               
                               
                               
                               lower.tail = FALSE) - x)
    
    
    
  }
  
  
  
  if (type=='mean'){ 
    
   if (min(c(mx))< 0){
       message("Estimated parameters produce an undefined mean.")
   }
   
    return(c(mx))

  }
  
  
  
  if (type== 'median'){
    
    
    
    return(c(px(.5)))
    
    
    
  }
  
  
  
  if (type == 'percentile'){
    
    
    
    return(c(px(p)))
    
    
    
  }
  
  
  
  if (type == 'all'){
    
    if(min(c(mx))< 0){
         message("Estimated parameters produce an undefined mean.")
     }
    
    return(data.frame(mean = mx, median = px(.5), percentile = px(p)))
    
    
    
  }
  
  
  
  else{
    
    
    
    return('invalid type')
    
    
    
  }
  
  
}
