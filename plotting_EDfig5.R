### model the population dynamics of an O12/O12-2 switching population
### Markus Arnoldini, October 2019

############### LOAD PACKAGES ###############
library(RColorBrewer)
library(ggplot2)
library(scales)
library(gridExtra)

############### DEFINE COLOR PALETTE ###############
mycols <- brewer.pal(2,name="Set2")

############ LOAD TIME DEPENDENT SIMULATION DATA ###############
det12    <- read_csv("det12.csv") 
det122   <- read_csv("det122.csv")

############ PLOT GROWTH CURVES, USE LINEAR Y AXIS FOR FINAL PRESENTATION ###############
#these are extended data figure 5D
p12 <- ggplot(det12,aes(x=time,y=O12,col='a')) +
  ggtitle("starting with O12 only") +
  geom_line() +
  geom_line(aes(x=time,y=O122,col='b')) +
  #scale_y_log10()+
  xlab("time (h)") +
  ylab("population size") +
  scale_color_manual(name="O antigen", values=c('a'=mycols[1],'b'=mycols[2]), labels = c("O12","O12-2")) +
  theme_bw()

p122 <- ggplot(data=det122,aes(x=time,y=O12,col='a')) +
  ggtitle("starting with O12-2 only") +
  geom_line() +
  geom_line(data=det122,aes(x=time,y=O122,col='b')) +
  #scale_y_log10()+
  xlab("time (h)") +
  ylab("population size") +
  scale_color_manual(name="O antigen", values=c('a'=mycols[1],'b'=mycols[2]), labels = c("O12","O12-2")) +
  theme_bw()

############ LOAD PARAMETER SCAN DATA ###############
scan12  <- read_csv("scan12start.csv") 
scan122 <- read_csv("scan122start.csv") 

############ PLOT PARAMETER SCANS ############### 
#these are extended data figure 5B
heatplot.122start <- ggplot(scan122,aes(x=s12_2,y=s2_12)) +
  theme_bw() +
  ggtitle("starting with O12-2 only") +
  xlab("switching rate from O12 to O12-2 (1/h)") +
  ylab("switching rate from O12-2 to O12 (1/h)") +
  geom_raster(aes(fill=fractionO12_2)) +
  scale_fill_gradientn(colors=c("yellow","blue"), name="fraction O12-2") +
  geom_contour(aes(z=fractionO12_2),col="white",linetype="dashed",size=0.3) 

heatplot.12start <- ggplot(scan12,aes(x=s12_2,y=s2_12)) +
  theme_bw() +
  ggtitle("starting with O12 only") +
  xlab("switching rate from O12 to O12-2 (1/h)") +
  ylab("switching rate from O12-2 to O12 (1/h)") +
  geom_raster(aes(fill=fractionO12_2)) +
  scale_fill_gradientn(colors=c("yellow","blue"), name="fraction O12-2") +
  geom_contour(aes(z=fractionO12_2),col="white",linetype="dashed",size=0.3) 
################


############ FIGURING OUT THE CORRECT SWITCHING RATES ###############
# PREDICTED FRACTIONS OF O12-2
# if starting with all O12-2
predfrac.122 <- 0.5
# if starting with all O12
predfrac.12 <- 0.1
# allowed variation in the predicted fraction
wiggle <- 0.05

############### GET DATA FOR SWITCHING RATES THAT FIT PREDICTED FRACTIONS FROM BOTH INITIAL CONDITIONS ###############
switching.122start <- scan122 %>% filter(fractionO12_2 > (predfrac.122-wiggle) & fractionO12_2 < (predfrac.122+wiggle))
switching.12start <- scan12 %>% filter(fractionO12_2 > (predfrac.12-wiggle) & fractionO12_2 < (predfrac.12+wiggle))

tempdata <- as.tibble(switching.122start) %>% mutate(initcon = "O12-2 only")
data     <- as.tibble(switching.12start) %>% mutate(initcon = "O12 only") 
data     <- rbind(tempdata,data)

############### MAKE LOCAL LINEAR REGRESSION TO DETERMINE INTERSECTION OF SWITCHING RATES ###############
dat12  <- switching.12start %>% filter(s2_12 <= 0.2 & s12_2 <= 0.2) %>% mutate(initcon="O12 only")
dat122 <- switching.122start %>% filter(s2_12 <= 0.2 & s12_2 <= 0.2) %>% mutate(initcon="O12-2 only")
reddat <- rbind(dat12,dat122)

model12  <- lm(dat12$s2_12 ~ dat12$s12_2)
model122 <- lm(dat122$s2_12 ~ dat122$s12_2)

#calculate where they intersect
icpt1 <- model12$coefficients[1]; icpt2 <- model122$coefficients[1]
slp1  <- model12$coefficients[2]; slp2  <- model122$coefficients[2]

xcoord <- (icpt2-icpt1)/(slp1-slp2) #O12 to O12-2 switching rate
ycoord <- slp1*xcoord + icpt1       #O12-2 to O12 switching rate

############### PLOT THE CORRECT SWITCHING RATES ONLY ###############
#these are extended data figure 5C
switchplot <- ggplot(data,aes(x=s12_2,y=s2_12,col=initcon)) + 
  theme_bw() +
  geom_rect(xmin=0,ymin=0,xmax=0.2,ymax=0.2,col=NA,fill='grey85',alpha=0.2) +
  geom_point(shape=1) +
  xlab("switching rate from O12 to O12-2 (1/h)") +
  ylab("switching rate from O12-2 to O12 (1/h)") +
  ggtitle("switching rates for observed population composition", subtitle = paste("target fractions ", 100*predfrac.122, "% (inital O12-2 only) and ", 100*predfrac.12, "% (initial O12 only)",sep="")) +
  scale_color_manual(name="initial composition", values=c('O12 only'=mycols[1],'O12-2 only'=mycols[2]), labels = c("O12 only","O12-2 only")) +
  theme(plot.title = element_text(size = 12))

zoomswitch <- ggplot(reddat,aes(x=s12_2,y=s2_12,col=initcon)) + 
  theme_bw() +
  coord_cartesian(xlim=c(0,0.2),ylim=c(0,0.2)) +
  geom_point(shape=1) +
  xlab("switching rate from O12 to O12-2 (1/h)") +
  ylab("switching rate from O12-2 to O12 (1/h)") +
  geom_abline(intercept=model12$coefficients[1],slope=model12$coefficients[2],linetype='dashed',alpha=0.5,col=mycols[1]) +
  geom_abline(intercept=model122$coefficients[1],slope=model122$coefficients[2],linetype='dashed',alpha=0.5,col=mycols[2]) +
  scale_color_manual(name="initial composition", values=c('O12 only'=mycols[1],'O12-2 only'=mycols[2]), labels = c("O12 only","O12-2 only"))

#grid.arrange(heatplot.122start,heatplot.12start,switchplot,zoomswitch,p122,p12,nrow=3)

