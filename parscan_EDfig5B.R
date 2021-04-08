### model the population dynamics of an O12/O12-2 switching population
### Markus Arnoldini, October 2019

############### LOAD PACKAGES ###############
library(deSolve)
library(tidyverse)

###GENERAL DEFINITIONS
gr <- 2.05

# PARAMETER RANGE FOR SCAN
s12_2.range <- seq(0, 1, by=0.01)
s2_12.range <- seq(0, 1, by=0.01)
outtime     <- 16

############### DEFINE DIFFERENTIAL EQUATIONS THAT MAKE UP THE MODEL ###############
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

############### DEFINE TIME STEPS ###############
timesteps <- seq(0, outtime+0.2, by=0.2)  #time is in h


# 1. START WITH ALL O12-2

#RUN PARAMETER SCANS OVER BOTH SWITCHING RATES
for(i in 1:length(s12_2.range)){
  for(j in 1:length(s2_12.range)){
    
    ############### DEFINE PARAMETERS ###############
    parameters <- c(
      mu     = gr,               #growth rate of STm, regardless of O-antigen, in 1/h
      s12_2  = s12_2.range[i],    #switching rate from O12 to O12-2, in 1/h, from probability 0.6 to switch per division
      s2_12  = s2_12.range[j],    #switching rate from O12-2 to O12, in 1/h, from probability 0.9 to switch per division
      CC     = 1E9                #carrying capacity
    )
    
    ############### DEFINE INITIAL CONDITIONS ###############
    initcons <- c(
      O12  = 0,      #population size of O12 STm
      O122 = 1E4     #population size of O12-2 STm
    )
    
    ############### SIMULATE ###############
    thisout <- data.frame(ode(y=initcons, times=timesteps, func=eqs, parms=parameters))
    
    ############### SELECT DATA AT DESIRED TIME POINT ###############
    t     <- outtime
    O12   <- thisout$O12[thisout$time==t]
    O12_2 <- thisout$O122[thisout$time==t]
    s12_2 <- s12_2.range[i]
    s2_12 <- s2_12.range[j]
    
    newout <- data.frame(t, O12, O12_2, s12_2, s2_12)
    
    if(i==1 && j==1){out <- newout}
    else{out <- rbind(out, newout)}
  }
}

############### ADD O12-2 RATIO TO THE OUTPUT DATA FRAME ###############
out.122 <- out %>% mutate(fractionO12_2 = O12_2/(O12+O12_2))

############### SAVE ###############
write.csv(out.122, "scan122start.csv")


# 2. START WITH ALL O12

#RUN PARAMETER SCANS OVER BOTH SWITCHING RATES
for(i in 1:length(s12_2.range)){
  for(j in 1:length(s2_12.range)){
    
    ############### DEFINE PARAMETERS ###############
    parameters <- c(
      mu     = gr,               #growth rate of STm, regardless of O-antigen, in 1/h
      s12_2  = s12_2.range[i],    #switching rate from O12 to O12-2, in 1/h, from probability 0.6 to switch per division
      s2_12  = s2_12.range[j],    #switching rate from O12-2 to O12, in 1/h, from probability 0.9 to switch per division
      CC     = 1E9                #carrying capacity
    )
    
    ############### DEFINE INITIAL CONDITIONS ###############
    initcons <- c(
      O12  = 1E4,  #population size of O12 STm
      O122 = 0     #population size of O12-2 STm
    )
    
    ############### SIMULATE ###############
    thisout <- data.frame(ode(y=initcons, times=timesteps, func=eqs, parms=parameters))
    
    ############### SELECT DATA AT DESIRED TIME POINT ###############
    t     <- outtime
    O12   <- thisout$O12[thisout$time==t]
    O12_2 <- thisout$O122[thisout$time==t]
    s12_2 <- s12_2.range[i]
    s2_12 <- s2_12.range[j]
    
    newout <- data.frame(t, O12, O12_2, s12_2, s2_12)
    
    if(i==1 && j==1){out <- newout}
    else{out <- rbind(out, newout)}
  }
}

############### ADD O12-2 RATIO TO THE OUTPUT DATA FRAME ###############
out.12 <- out %>% mutate(fractionO12_2 = O12_2/(O12+O12_2))

############### SAVE ###############
write.csv(out.12,"scan12start.csv")
