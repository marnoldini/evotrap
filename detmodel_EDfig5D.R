### model the population dynamics of an O12/O12-2 switching population
### Markus Arnoldini, October 2019

############### LOAD PACKAGES ###############
library(deSolve)
library(tidyverse)

############### DEFINE TIME STEPS ###############
timesteps <- seq(0, 24, by=0.2)  #time is in h

############### DEFINE DIFFERENTIAL EQUATIONS ###############
eqs <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    #O12 population size, logistic growth
    dO12  <- O12*(mu - mu*(O12+O122)/CC) - O12*(s12_2 - s12_2*(O12+O122)/CC) + O122*(s2_12 - s2_12*(O12+O122)/CC)
    #O12-2 population size, logistic growth
    dO122 <- O122*(mu - mu*(O12+O122)/CC) + O12*(s12_2 - s12_2*(O12+O122)/CC) - O122*(s2_12 - s2_12*(O12+O122)/CC)
    #return as list  
    list(c(dO12, dO122))
  })
}
############### DEFINE PARAMETERS ###############
parameters <- c(
  mu     = 2.05,      #growth rate of STm, regardless of O-antigen, in 1/h
  s12_2  = 0.0365,    #switching rate from O12 to O12-2, in 1/h, from probability 0.6 to switch per division
  s2_12  = 0.144,     #switching rate from O12-2 to O12, in 1/h, from probability 0.9 to switch per division
  CC     = 1E9        #carrying capacity
)

# 1. START WITH ALL O12

############### DEFINE INITIAL CONDITIONS ###############
initcons <- c(
  O12  = 1E4,      #population size of O12 STm, 
  O122 = 0         #population size of O12-2 STm
)
############### RUN ###############
out <- data.frame(ode(y=initcons, times=timesteps, func=eqs, parms=parameters))
############### SAVE ###############
write.csv(out,"det12.csv")

# 2. START WITH ALL O12-2

############### DEFINE INITIAL CONDITIONS ###############
initcons <- c(
  O12  = 0,        #population size of O12 STm, 
  O122 = 1E4       #population size of O12-2 STm
)
############### RUN ###############
out <- data.frame(ode(y=initcons, times=timesteps, func=eqs, parms=parameters))
############### SAVE ###############
write.csv(out,"det122.csv")
