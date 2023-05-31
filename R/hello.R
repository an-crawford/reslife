# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#source("gamma.R")

#We can combine all the functions in this file. For example:
#fsr = flexsurvreg(Surv(recyrs, censrec)~ 1, data = bc, dist = 'gamma')
#gamma.rl(fsr, 10, .99, 'all')

residlife <- function(data, x, p=.5, type = 'all', dist = ' ') {
  if (dist ==' '){
    return ('Please Specify Distribution!')
  }
  if (dist == 'gamma'){
    return(gamma.rl(data, x, p, type))
  }
  #add in other if statements
  else{
    return('Invalid Distribution!')
  }
  print("Hello, world!")
  print("Hello, world")



}
